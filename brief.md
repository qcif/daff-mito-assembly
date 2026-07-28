# Pipeline Specification: Organelle Genome Assembly + Barcode Recovery for Biosecurity Identification

**Status:** Draft v0.3 · **Format:** Nextflow DSL2 · **Downstream consumer:** Taxodactyl

---

# Hard constraints

These are client requirements and are not open to revision.

## 1. Purpose

Assemble organelle genomes (mitochondria, chloroplast) *de novo* from Oxford
Nanopore (ONT) sequencing of unknown / intercepted samples, then extract
taxonomic barcode loci from the assembled organelles and emit marker-tagged
barcode sequences plus metadata for taxonomic assignment by Taxodactyl (a
workflow we already built). The pipeline operates without prior knowledge of
the sample *species*, but assumes the submitter declares the target *kingdom*.

## 2. Scope

**In scope**
- ONT long reads (R10.4.1 / Dorado SUP basecalling assumed).
- Skim sequencing (~2× nuclear-equivalent). Organelles remain assemblable at
  this depth because they are high-copy (typically 30–150× organelle coverage
  from 2× nuclear).
- Kingdom-gated organelle assembly: plant (mt + cp), animal (mt).
- *De novo*, species-agnostic assembly within the declared kingdom (no
  same-species reference required).
- Output: assembled + annotated organelle contigs + marker-tagged barcode
  FASTA + per-sequence metadata.

**Boundary**
- Pipeline ends at **barcode extraction** from the assembled organelle(s).
  Taxonomic assignment, database search, candidate scoring and phylogenetics
  are Taxodactyl's responsibility.

**Out of scope**
- Fungi. Animal + plant only.
- Nuclear ribosomal (rDNA / ITS / 18S / 28S) recovery.
- Bacteria (chromosomes/plasmids).
- Cross-kingdom recovery from a single submission. Submissions without a
  declared kingdom are rejected at intake.
- Taxonomic assignment itself.
- Retention of non-recruited reads.


---

# Flexible constraints

All following constraints have been inferred from the above, and are open to
scrutiny and revision to achieve the intended goal.


## 3. Design Principles

1. **Assembly target is a hard pre-assembly gate.** The submitter declares
   `assembly_target` at intake — one of `animal_mt`, `plant_pt`, `plant_mt`.
   This selects the organelle reference index, coverage limits, Flye
   genome-size hint, BLAST DB, and genetic-code tables for downstream
   validation. Submissions with no declared target are rejected upstream and
   never reach this pipeline. Each sample-sheet row targets exactly one
   organelle (see [spec §1a](spec/01-pipeline-flow.md#1a-engineering-constraints));
   a plant submission requiring both plastid and mitogenome assembly submits
   two rows sharing the same reads. The declaration also identifies *what is
   to be characterised*, which licenses treating the host matrix as
   background: this narrows the search space and prevents a high-coverage
   host fraction from swamping the assembly graph.
2. **Host exclusion is positive recruitment, not negative depletion.** The
   kingdom gate is implemented by recruiting reads *toward* the target's
   organelle reference panel, not by depleting reads *against* a host
   reference. Positive selection requires no host reference (the submitter
   may not know the precise host) and cannot discard a genuine target read
   for incidentally resembling the host. Because conserved loci recruit
   fuzzily across taxon boundaries, recruitment is coarse enrichment only —
   precise target/background separation happens at post-assembly per-contig
   binning, with coverage and ORF validation applied regardless.
3. **Non-recruited reads are discarded.** Off-kingdom / off-target material
   is not retained by the pipeline. If a misdeclared kingdom is suspected,
   resubmit with corrected metadata.
4. **Whole-organelle assembly, then targeted barcode extraction.** Two
   logically distinct stages: (a) assemble the organelle genome(s) for the
   declared kingdom, (b) extract loci from the assembly using a config-driven
   panel. Keeping these decoupled lets the barcode panel evolve with Taxodactyl
   without re-running assembly, and lets the assembled organelle itself be a
   useful output.
5. **Samples are single-target; minor contamination tolerated.** Submissions
   are expected to contain one target organism per declared kingdom. Minor
   background contamination (a co-extracted host fragment, incidental
   environmental DNA) may occur and must not corrupt the target assembly;
   the assembler and post-assembly binning must therefore tolerate low-level
   non-target material. Emission is the dominant kingdom-matched organelle
   assembly by coverage; any secondary same-kingdom assemblies are recorded
   in diagnostics but not emitted.
6. **Full organelle annotation is preferred where cheap.** Beyond the
   Taxodactyl-driven barcode panel, run a proper organelle annotator
   (MITOS2 for animal mt, GeSeq or equivalent for plant cp/mt) if the
   runtime and packaging cost is modest. The annotated organelle is a
   useful standalone output; the barcode panel remains the pipeline's
   contractual deliverable to Taxodactyl.
7. **Locus panel is config-driven from Taxodactyl's accepted-loci list.** The
   extraction panel must not hardcode loci; it is parsed from (or generated to
   match) Taxodactyl's mutable accepted-loci set so the two never drift.
8. **Negative clarity.** A degraded sample yielding no organelle assembly, or
   an assembly with no extractable locus, produces an explicit no-recovery
   result with diagnostics — never indistinguishable from a confident
   negative.

## 4. Inputs

Pipeline entry point is a **sample sheet** (`samples.csv`) plus a
**data-directory root** that FASTQ paths are resolved against. Each row is
one sample. Multiple samples per invocation are supported and processed in
parallel. Per-sample outputs are keyed by `sample_id`. See `plan.md` §0 for
the sheet schema and validation rules.

| Input | Description |
|-------|-------------|
| `--samplesheet` | Path to `samples.csv` — one row per (sample, target) pair. Required columns: `sample_id`, `assembly_target` (one of `animal_mt`, `plant_pt`, `plant_mt`), `reads` (one or more **relative** FASTQ paths, pipe-delimited). Per [spec §1a](spec/01-pipeline-flow.md#1a-engineering-constraints), each row assembles exactly one organelle — a plant sample requiring both plastid and mitogenome submits two rows sharing the same reads. Optional columns: `sample_info`, `sample_type`, `sample_receipt_date`, `storage_location` — carried through to per-sample `metadata.json` and `report.html`; not used by pipeline logic. |
| `--data-dir` | Root directory used to resolve the relative `reads` paths in the sample sheet. Absolute paths in the sheet are rejected. |
| `locus_panel` | Value channel parsed from Taxodactyl accepted-loci config (shared across samples) |
| `organelle_refs` | Bundle root containing per-target organelle reference indices (`animal_mt.mmi`, `plant_pt.mmi`, `plant_mt.mmi`) — configurable, versioned, shared across samples |

## 5. Outputs

Outputs are emitted **per sample** under `outdir/<sample_id>/` (keyed by
`samples.csv` `sample_id`), plus a top-level run manifest and cross-sample
report. See `plan.md` §0 for the layout.

| Output (per sample) | Description |
|--------|-------------|
| `organelle_assembly` (FASTA + GFF/GenBank) | Polished organelle contig(s) for the target taxon, with full annotation (see §3.6) |
| `barcodes` (FASTA) | Marker-tagged barcode sequences extracted from the assembly |
| `metadata` (JSON) | Per-sequence: sample ID, contig ID, coverage, marker type, declared `assembly_target`, extraction coordinates, ORF/validation status, tool + DB versions |
| `report.html` | Per-sample workflow report: QC audit trail, recruitment statistics, secondary same-target contigs (if any), assembly + annotation summary, extracted barcode summary, and **explicit no-assembly and no-barcode-recovered results** distinct from empty/null output |

| Output (run-level) | Description |
|--------|-------------|
| `run_manifest.json` | Samplesheet snapshot, reference bundle version, pipeline commit — provenance for the whole invocation |
| `run-report.html` | Cross-sample overview: sample count, per-sample success / no-recovery breakdown, links into each per-sample `report.html` |

## 6. Module Layout

| Module | Tool(s) | Function / key decision |
|--------|---------|-------------------------|
| `QC_FILTER` | chopper, Filtlong, NanoPlot | Length/quality filter (chopper) + identity-weighted top-quality selection (Filtlong) + pre/post read stats (NanoPlot). Minimal filtering; adapter trimming is done upstream by Dorado at basecalling. |
| `RECRUIT` | minimap2 (or equiv.) | Positive recruitment against the target-specific organelle reference (`${assembly_target}.mmi`) — §3.2. Emits recruited reads → `COVERAGE_GATE`; non-recruited reads discarded (§3.3). Coarse enrichment; precision deferred to `BIN_TARGET`. |
| `COVERAGE_GATE` | seqkit + seqtk | Estimate coverage from recruited bases against the target's nominal organelle size (fixed per `assembly_target`). Subsample to configured MAX if over, **soft-fail** the sample if under configured MIN. Soft-fail is emitted as data (a status marker), not a pipeline error — sibling samples in the same run are unaffected. See `plan.md` §2.1 for limits and semantics. |
| `ASSEMBLE` | metaFlye (`--meta --nano-hq`) on gated reads (see `plan.md`) | Kingdom-routed organelle assembly on gated reads. Skipped for soft-failed samples. |
| `POLISH` | medaka | **Opt-in** (`--polish` boolean flag, default **off**). Homopolymer indel cleanup; enable when downstream ORF validation is failing on marginal-quality assemblies. Skipped by default because Dorado SUP R10.4.1 reads are already high-Q and medaka adds meaningful runtime. See `plan.md` §2 stage 8. |
| `BIN_TARGET` | coverage + ORF + reference comparison | Per-contig binning: separate true target organelle contigs from recruitment carry-over and low-level contamination using coverage spike, ORF integrity, and reference identity. Select dominant target by coverage; record any secondaries in diagnostics. |
| `ANNOTATE` | MITOS2 (animal mt) / GeSeq (plant cp+mt) | Full organelle annotation (§3.6). Feeds diagnostic visualisation and the barcode panel extraction. |
| `EXTRACT` | miniprot | Locate-and-extract each locus in `locus_panel` (protein-to-genome, works on divergent taxa with no same-species reference). Validate length, ORF, kingdom-appropriate genetic code, internal stops. |
| `COLLATE` | — | Per-sample aggregation: organelle FASTA + annotation + barcode FASTA + metadata JSON + `report.html`. |
| `RUN_REPORT` | — | Cross-sample overview: emit run-level `run-report.html` + `run_manifest.json` after all per-sample bundles complete. |

## 7. Validation & Safeguards

- **numt rejection:** nuclear mitochondrial pseudogenes are filtered via ORF /
  internal-stop checks (correct genetic-code table per marker) combined with
  the coverage filter (numts are low-copy relative to true organelle). Less
  of a concern post-organelle-assembly than in a whole-sample assembly, but
  retained as a safeguard.
- **Contamination tolerance:** assembler and BIN_TARGET must tolerate
  low-level non-target reads (see §3.5). Secondary same-kingdom assemblies,
  if any, are recorded in diagnostics; emission is target-only.
- **Coverage gate:** post-recruitment coverage estimation with per-kingdom
  MIN/MAX limits. Over-covered samples are subsampled; under-covered samples
  are soft-failed and skip assembly. Soft-fails are per-sample data, not
  pipeline errors — a failed sample does not block sibling samples in the
  same run. See `plan.md` §2.1.
- **Negative clarity:** see §3.8. Three distinct negative states — low
  coverage (gate soft-fail), no organelle assembled, and organelle assembled
  but no locus extractable — are individually reported.

## 8. Open Decisions

1. **Assembler choice.** Provisionally resolved to metaFlye on recruited
   reads (see `plan.md` §2 and §8.1). Benchmark against a kingdom-specific
   organelle assembler (GetOrganelle, Oatk) remains scheduled; the decision
   can be revised if the benchmark surfaces a clear winner.
2. **Read-level fallback arm** for samples too degraded to assemble. Kingdom
   gate already narrows the search space; a minimap2-against-reference fallback
   could rescue marginal samples. *Decision pending: include in v1, or add
   post-validation.*
3. **Polish depth** — whether medaka alone is sufficient or homopolymer-specific
   post-processing is warranted for marginal coverage.

---

# Resources

- Anna's mtDNA workflow: https://github.com/AnnaSyme/organelle-assembly
- CLAW Snakemake chloroplast WF: https://github.com/aaronphillips7493/CLAW
   - Perhaps steal [Flye/Minimap2 params](https://github.com/aaronphillips7493/CLAW/blob/main/config.yml) for Ont input data
