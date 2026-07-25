# Development Plan: Organelle Genome Assembly + Barcode Recovery

**Status:** Draft v0.2 · **Companion:** [brief.md](brief.md) · **Format:** Nextflow DSL2

This document expands [brief.md](brief.md) §6 into a concrete stage-by-stage
pipeline plan and lays out the development phases. Where this plan and the
brief disagree, this plan reflects the more recent decisions and the brief
should be reconciled (see §6).

---

## 0. Input: sample sheet

Pipeline is invoked with **two required parameters**:

- `--samplesheet <path>` — path to `samples.csv` (see schema below).
- `--data-dir <path>` — root directory that all `reads` paths in the sheet are resolved against.

One row per sample; multiple samples per invocation are supported and run in
parallel. All downstream stages fan out per-sample; outputs are keyed by
`sample_id`.

**Schema:**

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Unique per row. Used as the directory name for per-sample outputs. Must match `[A-Za-z0-9_.-]+` — no whitespace, no path separators. |
| `kingdom` | yes | `plant` \| `animal`. Rows with any other value (including empty) are rejected at parse time ([brief.md §1](brief.md)). |
| `reads` | yes | One or more **relative paths** to ONT FASTQ files, **pipe-delimited**. Paths are resolved against `--data-dir` (e.g. sheet says `run17/A.fastq.gz\|run17/A_pass2.fastq.gz`, `--data-dir /data`, → `/data/run17/A.fastq.gz` etc.). Absolute paths in the sheet are rejected. Multiple files are concatenated for that sample before `NANOPLOT_RAW`. `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz` all accepted. |
| `sample_info` | no | Free-text sample information. Carried through to per-sample `metadata.json` and rendered in `report.html`. Ignored by pipeline logic. |
| `sample_type` | no | Sample type descriptor (submitter-defined vocabulary). Carried through to `metadata.json` and `report.html`. Ignored by pipeline logic. |
| `sample_receipt_date` | no | Date the sample was received. Accepted format: ISO 8601 `YYYY-MM-DD`. Parse-warned (not failed) on other formats; the raw string is passed through. |
| `storage_location` | no | Physical storage location of the sample (submitter-defined). Carried through to `metadata.json` and `report.html`. Ignored by pipeline logic. |

**Example:**

Invocation: `nextflow run main.nf --samplesheet samples.csv --data-dir /data/run17`

```csv
sample_id,kingdom,reads,sample_info,sample_type,sample_receipt_date,storage_location
INT-2026-0007,plant,INT-2026-0007.fastq.gz,dried leaf fragment,leaf,2026-07-14,Freezer B / shelf 3 / box 12
INT-2026-0008,animal,INT-2026-0008_a.fastq.gz|INT-2026-0008_b.fastq.gz,merged two flowcells,,whole specimen,2026-07-18,
```

Resolved: `/data/run17/INT-2026-0007.fastq.gz`,
`/data/run17/INT-2026-0008_a.fastq.gz`, `/data/run17/INT-2026-0008_b.fastq.gz`.

**Validation (stage 0):**

- Header row: required columns (`sample_id`, `kingdom`, `reads`) must be present; optional columns (`sample_info`, `sample_type`, `sample_receipt_date`, `storage_location`) are recognised if present but not required. Order-flexible; case-sensitive. Unknown columns fail parse (guards against typos silently dropping fields).
- `sample_id` uniqueness enforced across the sheet; duplicates fail parse.
- `reads` paths must be relative — any leading `/` fails the row with an explicit "absolute path not allowed; use `--data-dir` root" error.
- All paths in `reads` are resolved against `--data-dir` and existence-checked before any sample starts. Missing files fail the whole run — no partial execution.
- `kingdom` is normalised to lowercase and matched against the enum. Anything else → parse error naming the offending row.
- `sample_receipt_date`, if present, is parsed as ISO 8601 (`YYYY-MM-DD`); malformed values emit a warning and pass through as the raw string.
- Optional columns are individually optional (empty cell is fine); an omitted column in the header simply means every sample has an implicit empty value for it.
- `--data-dir` must exist and be readable; failure to resolve is a fatal pre-flight error.

**Fan-out:** Nextflow `splitCsv(header: true)` on the samplesheet produces
one `tuple(meta, reads)` per row; `meta` carries `sample_id`, `kingdom`, and
any submitter-supplied optional fields (`sample_info`,
`sample_type`, `sample_receipt_date`, `storage_location`). All stages take
`meta` as a first-class channel key so per-sample context is threaded through
to `COLLATE`, where it is written into `metadata.json` and surfaced in
`report.html`.

**Per-sample output layout:**

```
outdir/
├── INT-2026-0007/
│   ├── organelle_assembly.fasta
│   ├── organelle_annotation.gff
│   ├── barcodes.fasta
│   ├── metadata.json
│   └── report.html
├── INT-2026-0008/
│   └── ...
├── run_manifest.json     # samplesheet snapshot + reference bundle version + pipeline commit
└── run-report.html  # cross-sample summary: N samples, success/no-recovery breakdown
```

## 1. Pipeline flow

FastP - long - alternative for QC? Ont suite alternatives?

```
samples.csv
   │
   ▼
[0] PARSE_SAMPLESHEET ─────────► one channel item per sample: (meta, reads)
   │                              meta = {sample_id, kingdom, sample_info}
   ▼    per-sample fan-out from here
[1] NANOPLOT_RAW ──────────────► raw read stats (HTML/PNG)
   │                              (adapters already trimmed by Dorado)
   ▼
[2] CHOPPER ──────────────────► quality-/length-filtered reads
   │
   ▼
[3] FILTLONG ─────────────────► top-quality reads (identity-weighted)
   │
   ▼
[4] NANOPLOT_CLEAN ────────────► post-clean read stats
   │
   ▼
[5] RECRUIT (minimap2) ────────► recruited reads (non-recruited discarded)
   │
   ▼
[6] COVERAGE_GATE ────────────► estimate coverage; subsample if > MAX;
   │                              fail sample (soft) if < MIN
   ▼
[7] METAFLYE (--meta --nano-hq)  draft organelle assembly + assembly graph
   │
   ▼
[8] MEDAKA (opt-in, --polish) ► polished assembly (skipped by default)
   │
   ▼
[9] BANDAGE_NG ───────────────► assembly graph image (diagnostic)
   │
   ▼
[10] BIN_TARGET (coverage + ref) target contig(s) + secondaries (diagnostic)
   │
   ▼
[11] BLAST_VALIDATE ──────────► per-contig identity to kingdom organelle refs
   │
   ▼
[12] ANNOTATE ────────────────► full organelle annotation (GFF/GenBank)
   │                              MITOS2 (animal) / TBD (plant, see §8)
   ▼
[13] MINIPROT_EXTRACT ────────► barcode coordinates + FASTA + ORF validation
   │
   ▼
[14] ORGANELLE_MAP ──────────────────► annotated organelle diagram (PNG/SVG)
   │
   ▼
[15] COLLATE ─────────────────► per-sample bundle under outdir/<sample_id>/
                                  (organelle FASTA + annotation, barcodes
                                   FASTA, metadata JSON, per-sample report)
   │
   ▼
[16] RUN_REPORT ──────────────► cross-sample summary (run-report.html)
                                  + run_manifest.json (samplesheet snapshot,
                                    ref bundle version, pipeline commit)
```

## 1a. Engineering constraints

Cross-cutting rules that every stage must satisfy. These are not
tool-specific; they apply uniformly to the whole workflow.

- **Every process runs in a container.** No stage may rely on tools
  installed on the host or in a conda env. Each `modules/local/*.nf`
  declares a `container` directive (or is resolved via a central
  `conf/containers.config` keyed by process name); the CI test run must
  fail if any process is missing one. Rationale: reproducibility across
  laptop → HPC → cloud, and clean auditability of tool provenance for
  biosecurity reporting. Applies equally to `COLLATE` and `RUN_REPORT`
  (Python + Jinja report renderer) — the report module ships as its own
  small image, not run against the host Python.
- **Pin image tags, never `latest`.** Version pinning is what makes
  `versions.yml` in the report trustworthy.
- **Prefer biocontainers or nf-core-published images.** Fall back to a
  minimal per-tool Dockerfile in `containers/<tool>/` only when no
  suitable published image exists.
- **File inputs must be stageable, not raw path strings.** On a remote
  executor (AWS Batch, Azure, Kubernetes, HPC scratch) Nextflow only
  copies a file to the compute node if it can see a Nextflow `Path`
  object referring to it. Two patterns produce one:
  1. **Preferred — declare as an `input: path` channel.** Materialise
     the param into a value channel at the top of `main.nf`
     (`ch_blastdb = Channel.value(file(params.blastdb))`) and pass it as
     `input: path blastdb` on every consuming process. This makes the
     dependency explicit in the process signature, keeps `script:`
     blocks free of `params.*` references, and lets a process be tested
     in isolation with a fixture path.
  2. **Acceptable — inline `${file(params.x)}` in the `script:` block.**
     `file()` returns a stageable `Path`, so this also works on remote
     executors:
     ```groovy
     blastdbcmd -db ${file(params.blastdb)} -entry_batch ${entry_batch} ...
     ```
     Use when the file is a stable, workflow-wide reference (e.g. a
     BLAST DB) and threading it through a channel would add noise
     without adding testability.

  What **does not** work: bare `${params.blastdb}` in a script block —
  that interpolates to a plain string, and the executor has nothing to
  stage. Also never write to paths outside the process work directory,
  and don't read from `$projectDir`/`$baseDir` at runtime for data.

  Rule of thumb: per-sample data and any file that fans out always go
  through explicit channels; workflow-wide reference bundles listed in
  [§4.4](#44-consolidated-build-script) may use either pattern. If in
  doubt, prefer pattern 1 — it's harder to misuse.

## 2. Stage detail

| # | Stage | Tool | Inputs | Outputs | Ref data | Key params / notes |
|---|-------|------|--------|---------|----------|--------------------|
| 1 | `NANOPLOT_RAW` | NanoPlot | raw FASTQ | read length/quality plots, summary TSV | — | Baseline metrics; feeds the report. Adapter/barcode trimming is performed upstream by Dorado at basecalling time, so no separate trim stage is needed. |
| 2 | `CHOPPER` | chopper | raw FASTQ | length-/quality-filtered FASTQ | — | Min length and min mean Q from config. Minimal filtering (preserve target reads). ONT-native. |
| 3 | `FILTLONG` | Filtlong | chopper-filtered FASTQ | identity-weighted top-quality FASTQ | — | Complements chopper's threshold cut with a percentile-based, identity-weighted selection (Filtlong scores reads by window quality and length). Prevents low-identity long reads from swamping METAFLYE at skim depth. Config: `--min_length`, `--keep_percent` (default 95). |
| 4 | `NANOPLOT_CLEAN` | NanoPlot | filtlong output | post-clean read stats | — | Pre/post comparison goes into report. |
| 5 | `RECRUIT` | minimap2 + samtools | cleaned FASTQ, kingdom organelle reference panel | `recruited.fastq.gz` | [§4.1](#41-setup-task-source-reference-panel-from-getorganelle) | Positive recruitment ([brief.md §3.2](brief.md)). Command shape: `minimap2 -ax map-ont -t $T $panel $reads \| samtools view -b -F 4 -q 1 -@ $T` — extract mapped read IDs, then `seqtk subseq` back to FASTQ. Liberal threshold — coarse enrichment only. Non-recruited reads discarded ([brief.md §3.3](brief.md)). Pattern adapted from [CLAW](CLAW/Snakefile). **Optional stricter filter** (P1 benchmark, plant-cp candidate): PAF-based recruitment with `minimap2 -cx map-ont`, retaining reads where alignment length / read length ≥ 0.7 and alignment length ≥ 1 kb — the [ptGAUL](reference-material/ptgaul/ptGAUL.sh) filter. Trades a small recall hit for cleaner input to METAFLYE; worth evaluating on plant samples where plastid recruitment is high-signal. |
| 6 | `COVERAGE_GATE` | seqkit stats + seqtk sample | recruited FASTQ, kingdom coverage limits (§2.1) | `gated.fastq.gz` (subsampled if needed), `coverage.json`, sample-status marker | — | Sum recruited bases via `seqkit stats -T`; estimate coverage as `total_bases / nominal_organelle_size`. If `< MIN` → **soft-fail** the sample: emit no-recovery marker, skip downstream, sample proceeds to `COLLATE` with an explicit low-coverage `report.html` (§2.1.3). If `> MAX` → `seqtk sample -s $seed $recruited $frac` where `frac = MAX / estimated`. If in-band → passthrough. Cross-sample isolation: `errorStrategy 'ignore'` on this process; failure state is data (a marker file), not a Nextflow error, so other samples are unaffected. |
| 7 | `METAFLYE` | Flye `--meta` (`--nano-hq` \| `--nano-corr`) | gated FASTQ | `assembly.fasta`, `assembly_graph.gfa`, `assembly_info.txt` | — | Retain `assembly_info.txt` — per-contig coverage is load-bearing for `BIN_TARGET`. `--meta` retained to tolerate low-level contamination ([brief.md §3.5](brief.md)). **Read-mode flag depends on whether an error-correction stage precedes Flye** — see §9 item 5. Default assumption: raw Dorado SUP → `--nano-hq`. Pass an organelle-keyed `--genome-size` hint (see §3.4). For plant samples, pass `--asm-coverage $((MAX-20))` per §3.5 (ptGAUL heuristic). |
| 8 | `MEDAKA` | medaka_consensus | assembly + gated reads | polished FASTA (only when enabled) | — | **Opt-in**, gated by `--polish` boolean flag (default `false` — polish is skipped unless the flag is passed). Homopolymer indel cleanup; useful when downstream ORF validation is failing on marginal-quality assemblies. Model matched to basecaller (Dorado SUP R10.4.1). By default the raw METAFLYE assembly flows straight into `BANDAGE_NG` / `BIN_TARGET`; the sample's `metadata.json` and `report.html` record the `polished` flag so ORF outcomes are attributable. Rationale for default-off: Dorado SUP R10.4.1 reads are already ~Q20+, medaka adds meaningful runtime, and the extra cleanup is not always needed — pass `--polish` when it is. |
| 9 | `BANDAGE_NG` | BandageNG | `assembly_graph.gfa` | graph image (PNG/SVG) | — | Diagnostic only; helps human review of tangled assemblies / multi-organism cases. |
| 10 | `BIN_TARGET` | minimap2 + custom | polished FASTA, `assembly_graph.gfa`, `assembly_info.txt`, kingdom refs | target contig(s) FASTA, secondaries TSV, plastid isoforms (plants only) | [§4.1](#41-setup-task-source-reference-panel-from-getorganelle) | Per-contig binning by coverage spike ∩ reference identity ∩ ORF integrity. Selects dominant target by coverage; secondaries logged not emitted ([brief.md §3.5](brief.md)). **Plant cp branch:** apply the quadripartite canonicalisation algorithm from §3.6 to emit `path1.fasta` (primary) + `path2.fasta` (alternative isoform); record edge-count structural signal in report. Ported from [ptGAUL](reference-material/ptgaul/combine_gfa.py). |
| 11 | `BLAST_VALIDATE` | blastn | target contigs, kingdom organelle reference DB | per-contig top hits + identity/coverage TSV | [§4.2](#42-refseq-organelle-dbs-for-blast_validate) | Sanity check that selected contigs match expected organelle type for the declared kingdom. Feeds metadata + report. |
| 12 | `ANNOTATE` | MITOS2 (animal mt) / **TBD** (plant cp+mt — see [§8](#8-remaining-open-questions)) | target contigs | annotation GFF + GenBank | — | Full organelle annotation ([brief.md §3.6](brief.md)). Feeds `ORGANELLE_MAP` and provides annotated organelle as a standalone output. Kingdom-dispatched: animal → MITOS2; plant → tool TBD. |
| 13 | `MINIPROT_EXTRACT` | miniprot | target contigs, locus panel (protein refs from Taxodactyl accepted-loci) | barcode FASTA, coordinates GFF, per-locus validation | [§4.3](#43-protein-panel-for-miniprot_extract) | Protein-to-genome; works on divergent taxa with no same-species reference. Validates length, ORF, kingdom-appropriate genetic code, internal stops ([brief.md §7](brief.md)). Runs independently of ANNOTATE — the barcode panel is the contractual output; ANNOTATE is supplementary. |
| 14 | `ORGANELLE_MAP` | **TBD** (see [§8](#8-remaining-open-questions)) | annotated GenBank from `ANNOTATE` | organelle map (PNG/SVG) | — | Diagnostic visualisation from the full annotation. Renders a circular/linear feature map for both plant (cp, mt) and animal (mt) targets. Tool choice deferred. |
| 15 | `COLLATE` | — | all upstream outputs (per sample) | `outdir/<sample_id>/{organelle_assembly.fasta, organelle_annotation.gff, barcodes.fasta, metadata.json, report.html}` | — | Per-sample aggregation; produces the Taxodactyl handoff bundle ([brief.md §5](brief.md)). Handles soft-failed samples: if the `COVERAGE_GATE` marker says fail, emit a minimal bundle with a low-coverage `report.html` and no assembly outputs. |
| 16 | `RUN_REPORT` | — | all `COLLATE` outputs + samplesheet | `outdir/run_manifest.json`, `outdir/run-report.html` | — | Cross-sample summary and manifest. Records samplesheet snapshot, reference bundle version (§4.4 `manifest.json`), pipeline commit, per-sample success/no-recovery/low-coverage breakdown. |

## 2.1 Coverage gate (§2 stage 6)

Post-recruitment, before assembly, estimate coverage from recruited bases and
either subsample or soft-fail. Purpose: (a) keep metaFlye out of the
over-coverage regime where the graph tangles, and (b) fail fast when there is
too little signal to recover an organelle, without letting a bad sample block
a batch.

### 2.1.1 Per-kingdom / per-organelle limits

Config-driven table (default values below; overridable via
`params.coverage_limits.<kingdom>.<organelle>`):

| Kingdom | Organelle | Nominal size (for coverage denominator) | MIN cov | MAX cov | Notes |
|---|---|---|---|---|---|
| plant | cp (primary) | 150 kb | 30× | 500× | Plants recruit cp + mt into one pool; cp is the dominant fraction (higher copy number, smaller genome). Gate on cp coverage — if cp passes, mt has enough to attempt too. |
| animal | mt | 17 kb | 30× | 300× | Tight organelle, high per-cell copy — MAX is generous because Flye handles small circles well up to a few hundred ×. |

Rationale for two-per-kingdom rather than one-per-kingdom: plant cp and animal
mt have very different nominal sizes and expected recruitment yields; a
single per-kingdom limit would over- or under-subsample one of them. Plant mt
does not get its own limit because it is not gated on directly — see the cp
note above.

### 2.1.2 Estimation formula

```
estimated_cov = total_recruited_bases / nominal_organelle_size
```

`total_recruited_bases` from `seqkit stats -T` on the RECRUIT output. No
mapping-based coverage estimate here — we don't have an assembly yet, and
counting bases against the nominal size is the standard pre-assembly proxy.

### 2.1.3 Decision matrix

| Condition | Action | Downstream effect |
|---|---|---|
| `estimated_cov < MIN` | **Soft-fail.** Write a `sample_status.json` marker with `status: "low_coverage"`, `estimated_cov`, `min_required`. Skip METAFLYE through ORGANELLE_MAP for this sample. | `COLLATE` sees the marker, emits a minimal per-sample bundle (`metadata.json` + `report.html` explaining the low-coverage failure). `RUN_REPORT` counts the sample under "low_coverage" in the run summary. |
| `MIN <= estimated_cov <= MAX` | Passthrough. Recruited FASTQ becomes `gated.fastq.gz` unchanged. | Normal flow through METAFLYE. |
| `estimated_cov > MAX` | Subsample. `seqtk sample -s $seed $recruited $frac` where `frac = MAX / estimated_cov`. Emit `coverage.json` recording pre- and post-subsampling depth and the seed. | Normal flow through METAFLYE on the subsampled FASTQ. |

Random seed default `--seed 42`; overridable per-run via `params.seqtk_seed`
for reproducibility of downsampled inputs.

### 2.1.4 Cross-sample failure isolation

**Requirement:** a soft-failed sample must not block any other sample in the
same invocation.

Implementation:

- `COVERAGE_GATE` never exits non-zero on a coverage decision — MIN/MAX
  outcomes are *data*, written to `sample_status.json`. Exit code 0 in all
  three branches.
- Nextflow process is configured `errorStrategy 'ignore'` **only** to guard
  against unexpected tool failures (e.g. seqkit or seqtk crash); the coverage
  decision itself does not rely on this.
- Downstream stages (METAFLYE onward) are conditionally executed on
  `sample_status.status == "ok"`. Nextflow filter on the channel:
  `gated_ch.branch { ok: it.status == "ok"; failed: true }` — the failed
  branch skips straight to `COLLATE`.
- `COLLATE` accepts both branches; on the failed branch it emits the minimal
  bundle described in 2.1.3.
- `RUN_REPORT` is the join point across all samples; it does not care whether
  each sample succeeded or soft-failed, only that a `metadata.json` exists.

This means: **a single low-coverage sample in a batch of 100 does not stop
the other 99**, and the failure is fully visible in both the per-sample and
run-level reports rather than being an opaque pipeline error.

## 2.2 Custom logic components

Most stages are thin wrappers around off-the-shelf bioinformatics tools
(`minimap2`, `Flye`, `medaka`, `blastn`, `miniprot`, etc.) invoked
directly from a `script:` block. This section catalogues the places
where custom logic — Python scripts, non-trivial parsing, or novel
decision code — is required. Each item is:

- an in-house script living under `bin/` (Nextflow auto-stages `bin/`
  onto every process `PATH`);
- distributed inside a **small dedicated container** (`containers/<name>/`)
  per [§1a](#1a-engineering-constraints) — no host-Python fallback;
- covered by unit tests independent of the Nextflow harness, so the
  logic is testable without spinning up the pipeline.

If you find yourself writing more than ~10 lines of `awk`/`bash` inside
a `.nf` script, promote it into `bin/` and add an entry here.

| # | Component | Container | Stage(s) | Purpose | Spec anchor |
|---|---|---|---|---|---|
| C1 | `parse_samplesheet.py` | `wf5/samplesheet:<tag>` (thin Python + `pandas`/`csv` + `pyyaml`) | [§2 stage 0](#2-stage-detail) | CSV validation beyond nf-schema's reach: `sample_id` uniqueness across the sheet, pipe-split of `reads`, absolute-path rejection with an actionable error, resolution against `--data-dir`, existence check of every resolved path, kingdom normalisation, ISO 8601 date warning-not-fail. Emits a normalised JSONL that Nextflow `splitCsv`-consumes into the per-sample channel. | [§0](#0-input-sample-sheet) |
| C2 | `coverage_gate.py` | `wf5/coverage-gate:<tag>` (Python stdlib only + `seqkit`/`seqtk` in `PATH`) | [§2 stage 6](#2-stage-detail) | Reads `seqkit stats -T` output, computes `estimated_cov = total_bases / nominal_organelle_size` from the kingdom-keyed limits table, executes the three-branch decision (soft-fail / passthrough / subsample), invokes `seqtk sample -s $seed` when subsampling, and writes `sample_status.json` + `coverage.json`. **Always exits 0** so a coverage decision is data, not a Nextflow error ([§2.1.4](#214-cross-sample-failure-isolation)). | [§2.1](#21-coverage-gate-2-stage-6) |
| C3 | `bin_target.py` | `wf5/bin-target:<tag>` (Python + `biopython` + `pysam`/`mappy` for GFA/BAM parsing) | [§2 stage 10](#2-stage-detail) | Per-contig classification by the intersection of (a) coverage spike vs. background, (b) minimap2 identity to the kingdom organelle reference panel, and (c) ORF integrity under the kingdom-appropriate genetic code. Selects the dominant target contig, records secondaries in a diagnostic TSV, and for animal samples runs the end-overlap circularity check flagged in [§3.3](#33-specific-issues-and-decisions). | [§3.3](#33-specific-issues-and-decisions), [§9 item 10](#9-fine-tuning-post-prototype-benchmarking) |
| C4 | `plastid_canonicalise.py` (ported from [ptGAUL `combine_gfa.py`](reference-material/ptgaul/combine_gfa.py)) | shares `wf5/bin-target:<tag>` (invoked by C3 on the plant-cp branch) | [§2 stage 10](#2-stage-detail), plant-cp branch only | Parses Flye's `assembly_graph.gfa`, classifies edges by length + depth (LSC = longest, IR = deepest, SSC = remainder for the 3-edge case), emits `path1.fasta` and `path2.fasta` for the two SSC-orientation isoforms. Handles the 1-edge (resolved circle passthrough) and N≠3 (diagnostic-only) cases explicitly. **Ported inline** because upstream is unmaintained. | [§3.6](#36-plastid-quadripartite-canonicalisation-ptgaul-derived) |
| C5 | `validate_barcodes.py` | `wf5/barcode-validate:<tag>` (Python + `biopython`) | [§2 stage 13](#2-stage-detail) | Post-processes miniprot output: per-locus length check, ORF check under the kingdom-appropriate NCBI genetic-code table, internal-stop check, protein-identity threshold. Emits `barcodes.fasta` (validated only) + a per-locus `validation.tsv` recording pass/fail reasons for panel loci that dropped out. **Animal-mt clade trial** ([§3.3](#33-specific-issues-and-decisions)): attempts NCBI tables 2 (vertebrate) and 5 (invertebrate), picks the one with a valid ORF, records the chosen table in metadata. | [§3.2](#32-fork-vs-dynamic-parameter--by-stage), [§3.3](#33-specific-issues-and-decisions) |
| C6 | `collate.py` | `wf5/report:<tag>` (shared with C7 and the report renderer — Python + Jinja2 per [§6a.5](#6a5-report-implementation-boilerplate)) | [§2 stage 15](#2-stage-detail) | Detects the [`COVERAGE_GATE`](#21-coverage-gate-2-stage-6) soft-fail marker and dispatches to either the full or minimal per-sample bundle, aggregates upstream outputs into the layout in [§0](#0-input-sample-sheet), emits `metadata.json` (sample meta + tool versions + reference-bundle version), and invokes the [report renderer](wf-report-boilerplate/report.py) to write `report.html`. | [§6a.5](#6a5-report-implementation-boilerplate) |
| C7 | `run_report.py` | shares `wf5/report:<tag>` with C6 | [§2 stage 16](#2-stage-detail) | Cross-sample join: reads every per-sample `metadata.json`, classifies each sample into `ok` / `no_recovery` / `low_coverage` / `hard_failure`, emits `run_manifest.json` (samplesheet snapshot + reference-bundle version from [§4.4 manifest](#44-consolidated-build-script) + pipeline commit + invocation timestamp), and renders the run-level HTML via the same Jinja machinery as C6. | [§6a.4](#6a4-run-level-run-reporthtml) |

**Non-workflow custom logic (pipeline-adjacent):**

- `scripts/build_refs.sh` (or `.py`) — the consolidated reference-bundle
  builder described in [§4.4](#44-consolidated-build-script). Not
  invoked by the workflow at run time; produces the versioned
  `refs/v<YYYY.MM>/` directory that `params.kingdom_refs` points at.
  Emits `manifest.json` with SHA256 digests, source URLs, release tags
  per input. Runs in its own container (`wf5/ref-build:<tag>`) with
  `minimap2` + `blast+` + `datasets`/`efetch` + `pandas`.

**Container consolidation:** C3 and C4 share one container because C4 is
invoked as a helper by C3. C6 and C7 share one container with the report
renderer because all three are Python + Jinja2 with the same static
assets. C1, C2, and C5 each stand alone — their tool surfaces don't
overlap enough to justify a fatter shared image. Net: **five bespoke
images** for custom code, plus per-tool biocontainers for the
off-the-shelf stages (see [§1a](#1a-engineering-constraints)).

**Deps in the image, code at runtime.** Each bespoke image is a lean
base (Python interpreter + `pip install -r requirements.txt` + any
system tools the component shells out to). **The Python source is not
baked in.** At runtime, the source under `scripts/` is delivered to
the container so a code edit does not require a rebuild:

- **Nextflow-native path (all executors):** helper scripts live under
  `bin/` — Nextflow auto-stages `bin/` onto every process's `PATH` in
  the work directory. Any script placed there is callable by name
  inside the container. This works on Docker, Azure, AWS Batch, K8s,
  and HPC identically.
- **Local dev shortcut (Docker only):** for iterating on a
  larger `scripts/` tree without touching `bin/`, add
  `containerOptions '-v ${projectDir}/scripts:/opt/wf5/scripts:ro'`
  on the process (guarded behind a profile so it's dev-only). Do
  not use this on remote executors — the host path won't exist on
  the compute node.

Net effect: `images.yml` rebuilds only when `requirements.txt` or the
`Dockerfile` changes ([§5b](#5b-ci)); day-to-day Python edits ship
via the next `nextflow run` without a Docker step.

**Test surface:** each `bin/*.py` gets a unit-test file under
`tests/unit/` covering its decision logic against synthetic fixtures.
Nextflow-level integration testing is separate (P0 `-profile test`
covers channel wiring; P1+ tests exercise real tools). The unit tests
run in the component's own container in CI, matching the runtime
environment exactly.

## 3. Organelle considerations: mitochondrion vs plastid

Plants carry both a mitogenome and a plastome; animals carry only a
mitogenome. The two are similar enough that a single assembly pass handles
both, but several downstream parameters must vary by organelle type. This
section enumerates what differs and where the workflow forks.

### 3.1 What differs between organelles

| Aspect | plant cp | plant mt | animal mt |
|---|---|---|---|
| Typical size | 120–160 kb | 200 kb – several Mb | 15–20 kb |
| Structure | quadripartite, large inverted repeat | multipartite, recombinant, multiple isoforms | simple circle, compact |
| Relative coverage | very high (>1000× possible) | medium | medium–high |
| Genetic code (ORF validation) | NCBI table 11 (bacterial/plastid) | NCBI table 1 (standard) | table 2 (vertebrate) / table 5 (invertebrate) — clade-dependent |
| Primary barcodes | rbcL, matK, trnH-psbA, ndhF, atpB | rarely used for barcoding | COX1, CYTB, 12S, 16S |
| Annotation tool | TBD ([§8](#8-remaining-open-questions)) | TBD ([§8](#8-remaining-open-questions)) | MITOS2 (or MitoFinder) |
| BLAST validation DB | RefSeq plastid | RefSeq plant mitochondrion | RefSeq metazoan mitochondrion |
| Assembly pitfall | inverted repeat collapses or branches ambiguously — canonical quadripartite form not guaranteed (canonicalised post-assembly per §3.6) | metaFlye emits multiple alternative isoforms; no canonical form | should close cleanly as one circle; verify circularity |

### 3.2 Fork vs dynamic parameter — by stage

| Stage | Approach | Rationale |
|---|---|---|
| `RECRUIT` | **dynamic**: kingdom-keyed reference panel bundles all relevant organelles (plant = cp + mt in one panel) | A single recruitment pass keeps the read pool unified for assembly. Plant kingdom panel must include both cp and mt references. |
| `METAFLYE` | **shared**: one assembly pass per sample | No assembly-time fork. metaFlye handles mixed-organelle input fine; contigs separate naturally by coverage and gene content. |
| `MEDAKA` | **shared, opt-in**: one polish pass when `--polish` is passed | Identical model regardless of organelle. Skipped by default; enable with `--polish`. |
| `BANDAGE_NG` | **shared**: one graph image | The full graph shows cp + mt branching for plants — diagnostically useful. |
| `BIN_TARGET` | **fork point**: classify each contig as cp / mt / off-target by reference identity + coverage tier | Natural fork point. Downstream stages branch per-contig from here. |
| `BLAST_VALIDATE` | **dynamic**: per-organelle BLAST DB selected by `BIN_TARGET` label | Validation must match the contig's classified organelle. |
| `MINIPROT_EXTRACT` | **dynamic**: per-contig genetic code table and locus panel subset selected by organelle label (and clade-trial for animal mt — see §3.3) | Genetic code differs by organelle; locus panel partitions naturally (rbcL → cp, COX1 → mt). |
| `ANNOTATE` | **dynamic**: MITOS2 for animal mt, TBD for plant cp+mt ([§8](#8-remaining-open-questions)) (per organelle label from `BIN_TARGET`) | Annotator is organelle-specific. |
| `ORGANELLE_MAP` | **shared**: consumes annotation output from `ANNOTATE` | Renderer is generic over organelle type. |
| `COLLATE` | **shared**: aggregate all per-organelle outputs into one bundle | Plant samples produce cp + mt barcode FASTAs concatenated into a single emission. |

### 3.3 Specific issues and decisions

- **Plant samples emit two organelles per target organism.** Animal samples
  emit one (mt); plant samples emit both cp and mt of the same target plant.
  Brief.md §3.5 wording covers this via "dominant kingdom-matched organelle
  assembly" — for plants both organelles of the dominant organism are in
  scope.

- **Plastid inverted repeat.** metaFlye may collapse the IR or emit alternative
  paths; canonical quadripartite form (LSC–IRb–SSC–IRa) is not guaranteed. For
  barcode extraction this does not matter — miniprot finds genes regardless of
  IR resolution. For ORGANELLE_MAP visualisation it does; a concrete canonicalisation
  algorithm derived from [ptGAUL](reference-material/ptgaul/ptGAUL.sh) is
  specified in §3.6.

- **Plant mitogenome multi-isoform.** Plant mt recombines into multiple
  alternative arrangements; metaFlye emits multiple contigs with shared
  regions. Acceptable for barcoding — gene-bearing contigs are what matter.
  Diagnostics should flag when multiple alternative mt contigs are present.

- **Animal mt genetic code is clade-dependent.** Vertebrate (table 2),
  invertebrate (table 5), echinoderm (table 9), etc. The submitter declares
  kingdom but not phylum. Options: (a) default to invertebrate (table 5) —
  covers most biosecurity intercepts; (b) require phylum at intake; (c) try
  table 2 and table 5 in EXTRACT, pick the one with valid ORF.
  **Recommend (c)** — automatic, unambiguous when an ORF is recoverable, and
  avoids burdening the submitter.

- **rDNA out of scope** ([brief.md §2](brief.md)). Nuclear rDNA (ITS, 18S,
  28S) recovery is dropped from this workflow. Extraction is limited to
  organelle-encoded loci.

- **Animal mt circularisation check.** Animal mt should close as a single
  ~16 kb circle. End-overlap circularity check added to `BIN_TARGET` is
  diagnostically valuable; flag in P3.

### 3.4 Flye `--genome-size` hint (per organelle)

Flye accepts a `--genome-size` hint (e.g. `--genome-size 135k`) that improves
assembly on small, high-coverage targets like organelles.
[CLAW](CLAW/config.yml) hard-codes 135 kb for chloroplasts; we generalise per
organelle label emitted by `BIN_TARGET` — except Flye runs *before*
`BIN_TARGET`, so the hint at the METAFLYE stage is set per-**kingdom** using
the largest expected organelle in that kingdom (so we don't under-hint a
plant mitogenome).

| Kingdom | Hint at METAFLYE (largest expected) | Rationale |
|---|---|---|
| plant | `2m` | Plant mt can reach 200 kb–several Mb; cp is small enough not to be prejudiced. |
| animal | `20k` | Animal mt is ~15–20 kb; tight hint helps recover the small circle. |

The hint is advisory to Flye's coverage estimator, not a length filter — over-
or under-hinting degrades but does not silently truncate. Re-evaluate the
plant hint in P2 if mt assembly fragments.

### 3.5 Flye `--asm-coverage` for plastid runs

For plant samples, when a plastid contig is being targeted at high coverage,
pass Flye `--asm-coverage <N>` where `N = coverage_max − 20` (COVERAGE_GATE's
plant MAX minus 20; see §2.1.1). Rationale (adopted from
[ptGAUL](reference-material/ptgaul/ptGAUL.sh)): when total coverage exceeds
this value, Flye uses only the longest reads for the initial contig-building
step, which improves contiguity on high-copy small circles like plastids.
Skip for animal mt — the target is too small for `--asm-coverage` to matter.

### 3.6 Plastid quadripartite canonicalisation (ptGAUL-derived)

Plastid genomes have a canonical quadripartite structure: **LSC** (large
single-copy, ~80–90 kb) – **IRa** (inverted repeat A, ~20–30 kb) – **SSC**
(small single-copy, ~15–20 kb) – **IRb** (reverse-complement of IRa). Because
the IR is present twice in the plastome, a Flye assembly graph typically
emits **three edges**: LSC, SSC, and a single IR contig (whose read coverage
is roughly double the SC edges).

The two SSC orientations are biologically real — plastids exist as a mixture
of both isoforms — so a "correct" plastid assembly has two valid linear
representations (`path1` and `path2`). This is implemented in
[ptGAUL's combine_gfa.py](reference-material/ptgaul/combine_gfa.py) and is
lifted directly into our `BIN_TARGET` plant-cp branch.

**Algorithm:**

1. Parse Flye's `assembly_graph.gfa`; extract per-edge sequence (from `S` lines) and per-edge depth (from the `S` line depth tag).
2. Count edges. Classify:
   - **3 edges** → canonical quadripartite. Proceed to canonicalise.
   - **1 edge** → fully-resolved circle. Passthrough as the assembly; flag "no IR resolution needed" in the report.
   - **anything else** → non-canonical. Emit as-is with a diagnostic flag; do not attempt automated canonicalisation. Human review via BANDAGE_NG output.
3. For the 3-edge case:
   - Longest edge by sequence length → **LSC**.
   - Deepest edge by read depth (top of `sort -k2 -n -r`) → **IR** (2× coverage betrays it).
   - Remaining edge → **SSC**.
4. Emit two canonical isoforms:
   - `path1.fasta`: `LSC + IR + SSC + reverse_complement(IR)`
   - `path2.fasta`: `LSC + IR + reverse_complement(SSC) + reverse_complement(IR)`
5. Record both isoforms in per-sample outputs. Downstream stages (BLAST_VALIDATE, ANNOTATE, ORGANELLE_MAP, MINIPROT_EXTRACT) run against `path1` as the primary; `path2` is emitted alongside and noted in `metadata.json` as the alternative isoform.

**Edge-count QC signal** is reported explicitly in the per-sample
`report.html` under the "Assembly quality assessment" section (§6a.2):
"3 edges → canonical", "1 edge → resolved circle", or "N edges → manual
review recommended". This gives the operator an immediate structural sanity
check without needing to open BandageNG.

Implementation: port `combine_gfa.py` (small, ~90-line Python) into the
pipeline as a stage-internal script; upstream is unmaintained and pinning it
directly is cleaner than a conda dep.

## 4. Reference data

Kingdom-keyed reference bundle, versioned and configurable via `params.kingdom_refs`:

| Kingdom | Organelle refs | Protein panel (for miniprot) |
|---------|---------------|------------------------------|
| plant | RefSeq plastid + plant mitogenome | rbcL, matK, ndhF, atpB, COX1, etc. |
| animal | RefSeq metazoan mitogenome | COX1, CYTB, 12S, 16S, etc. |

Locus panel content is **not authoritative here** — it is parsed at runtime
from Taxodactyl's accepted-loci config ([brief.md §3.7](brief.md)). The table
above shows representative loci only.

### 4.1 Setup task: source reference panel from GetOrganelle

Decision: bootstrap the kingdom organelle reference panel from
[GetOrganelleDB](https://github.com/Kinggerm/GetOrganelleDB) (the reference
library shipped with the GetOrganelle assembler). The libraries are already
kingdom-partitioned and curated for organelle recruitment. We reuse them as
recruitment references without using GetOrganelle as the assembler.

**Mapping from GetOrganelle libraries to our kingdom panels:**

| Our kingdom | GetOrganelle libraries | Notes |
|---|---|---|
| plant | `embplant_pt` + `embplant_mt` | Plastid + mitogenome. Use `other_pt` if non-embryophyte plant lineages need coverage. |
| animal | `animal_mt` | Mitogenome only. |

**Setup procedure (one-off, scripted, versioned):**

1. Fetch a pinned [GetOrganelleDB release](https://github.com/Kinggerm/GetOrganelleDB/releases) tarball; record the release tag in `params.kingdom_refs.version`.
2. Extract per-library FASTA files.
3. Concatenate the libraries listed above into one FASTA per kingdom.
4. Build minimap2 indices: `minimap2 -d kingdom_<plant|animal>.mmi kingdom_<...>.fa`.
5. Stage indices + source manifests under a versioned directory (e.g. `refs/v2026.06/`); point `params.kingdom_refs` at it.
6. Record source URLs, release tags, and SHA256 digests in a `manifest.json` alongside the indices. Emit the manifest version into per-run metadata output ([brief.md §5](brief.md)).

**Re-evaluate after P1.** If GetOrganelle library coverage proves insufficient on the P1 plant test data (older or non-model lineages underrepresented), fall back to a self-built RefSeq-derived panel using the canonical sources listed in §4.

### 4.2 RefSeq organelle DBs (for `BLAST_VALIDATE`)

Stage 11 needs kingdom-appropriate BLAST nucleotide DBs, separate from the
recruitment panel — recruitment is fuzzy and permissive; validation must be
against a clean, curated reference set to catch mis-binning.

- **Source:** NCBI RefSeq organelle FTP:
  - plastid: `https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*.genomic.fna.gz`
  - mitochondrion: `https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*.genomic.fna.gz`
- **Kingdom split:** RefSeq mitochondrion is not pre-split by kingdom. Split by NCBI taxonomy at build time using accession → taxid mapping (`efetch` or the `nucl_gb.accession2taxid` dump) and the taxonomy dump (`taxdump.tar.gz`). Emit two files: `refseq_mt_metazoa.fa`, `refseq_mt_viridiplantae.fa`. Plastid is plants-only by definition.
- **Version pinning:** record RefSeq release number (e.g. `RefSeq 227`) + download date + taxdump date in `manifest.json`.
- **Build:** `makeblastdb -in refseq_<kingdom>_<organelle>.fa -dbtype nucl -parse_seqids -out blastdb/refseq_<...>`.
- **Refresh cadence:** re-fetch quarterly; RefSeq updates every ~2 months.

### 4.3 Protein panel for `MINIPROT_EXTRACT`

miniprot needs a protein FASTA per locus (protein-to-genome alignment). The
locus *identity* comes from Taxodactyl's accepted-loci config at runtime
([brief.md §3.7](brief.md)); the protein *sequences* are what we fetch and
stage here.

- **Source of truth for locus list:** Taxodactyl accepted-loci config, referenced by `params.locus_panel` (a filesystem path). The file is staged onto every consuming process via a value channel and parsed inside the container ([§1a](#1a-engineering-constraints) channel rule); it is never read from `params.*` inside a `script:` block or preloaded into a Groovy map at workflow entry. Never hardcoded here.
- **Source of protein sequences:** [NCBI RefSeq protein](https://ftp.ncbi.nlm.nih.gov/refseq/release/) — query by gene symbol + kingdom taxid restriction (e.g. `rbcL AND txid33090[Organism]` for plants). Alternative: [UniProt](https://www.uniprot.org/) reviewed entries, filtered by taxonomy.
- **Per-locus curation:** pull one representative protein per major clade within the target kingdom (e.g. for COX1: one per invertebrate order + one vertebrate representative). miniprot tolerates divergence — 5–10 representatives per locus is enough; more slows alignment without improving sensitivity.
- **Version pinning:** record per-locus accessions + fetch date in `manifest.json` under `locus_panel_proteins`.
- **Build:** one FASTA per locus per kingdom: `proteins/<kingdom>/<locus>.faa`. Concat per kingdom for a single miniprot invocation, or pass individually depending on runtime characteristics measured in P3.
- **Refresh trigger:** re-fetch when Taxodactyl's accepted-loci config changes, or annually, whichever comes first.

### 4.4 Consolidated build script

All three reference-build tasks (§4.1–§4.3) live in a single versioned script
(`scripts/build_refs.sh` or equivalent) that emits `refs/v<YYYY.MM>/` with:

```
refs/v2026.06/
├── recruit/
│   ├── plant.mmi         # GetOrganelle embplant_pt + embplant_mt
│   └── animal.mmi        # GetOrganelle animal_mt
├── validate/
│   ├── refseq_pt.{nhr,nin,nsq,...}      # plant plastid BLAST DB
│   ├── refseq_mt_metazoa.{...}
│   └── refseq_mt_viridiplantae.{...}
├── proteins/
│   ├── plant/<locus>.faa
│   └── animal/<locus>.faa
└── manifest.json         # SHA256 + version + source URL for every input
```

`manifest.json` is the single source of truth for reference provenance; its
version string is emitted into per-run metadata ([brief.md §5](brief.md)) so
every result is traceable to the exact reference bundle used.

## 5. Test data

Per kingdom, one clean + one contaminated dataset. Skim depth
(~2× nuclear-equivalent) throughout — see [brief.md §2](brief.md).

| Kingdom | Clean target | Minor-contamination scenario |
|---------|-------------|------------------------------|
| plant | known reference plant (e.g. *Arabidopsis*) ONT skim | plant + trace insect contaminant |
| animal | known insect (e.g. *Drosophila*) ONT skim | insect + trace plant host DNA |

Acceptance criteria per dataset:
- Clean: full organelle assembled, all panel loci extractable with valid ORF,
  full annotation produced by `ANNOTATE`.
- Contamination: target dominant assembly selected; any low-coverage secondary
  contigs recorded in diagnostics; target assembly unaffected by trace
  non-target reads.

## 5a. Tests

Two testing surfaces, one per language, both run in CI on every push
inside the same containers the workflow uses at runtime
([§1a](#1a-engineering-constraints)):

- **Python (`bin/*.py`, `wf-report-boilerplate/`) — `unittest`.**
  One test file per module under `tests/unit/`, imported directly (no
  Nextflow harness). Each [custom-logic component](#22-custom-logic-components)
  (C1–C7) gets its own test module. Fixtures live under
  `tests/fixtures/` — synthetic GFA files, minimal FASTA snippets,
  hand-authored `seqkit stats` outputs — small enough to check in.
  `flake8` runs alongside the tests (per user CLAUDE.md).
- **Nextflow (`main.nf`, `modules/local/*.nf`) — [nf-test](https://www.nf-test.com/).**
  Per-module tests exercise each process's channel wiring, `input:` /
  `output:` shapes, and stub behaviour in isolation. A workflow-level
  test replays the P0 `-profile test` end-to-end run and asserts the
  per-sample output layout in [§0](#0-input-sample-sheet). Reference
  bundles are pointed at the same `tests/fixtures/` used by the Python
  unit tests.

**Coverage target: 100%.** Where a component is genuinely untestable in
CI — because it is variable (calls out to a network service, is
non-deterministic), expensive (needs a real ONT dataset or a full
BLAST DB), or too complicated to fixture (`METAFLYE`, `MEDAKA`,
`ANNOTATE`, `ORGANELLE_MAP`) — mock the tool call at the process boundary and
test only our wrapper logic. Record every mocked-out surface in the
per-component test module's docstring so it's visible what is *not*
being covered.

**Rule of thumb for what to mock:**

- Off-the-shelf tools invoked by pure-wrapper processes
  (`NANOPLOT_*`, `CHOPPER`, `FILTLONG`, `METAFLYE`, `MEDAKA`,
  `BANDAGE_NG`, `BLAST_VALIDATE`, `ANNOTATE`, `ORGANELLE_MAP`,
  `MINIPROT_EXTRACT`) — mock at nf-test level; assert we assemble the
  command line correctly and consume the output shape correctly. Do
  not try to actually run the tool.
- Custom-logic components (C1–C7 in [§2.2](#22-custom-logic-components))
  — real Python unit tests with real fixtures. These are our code and
  they must be covered.
- Network fetches (`scripts/build_refs.sh` calls to NCBI FTP) — mock
  the fetch, real code for parsing.

CI matrix reports coverage per component; a drop below 100% for
non-mocked lines fails the build.

## 5b. CI

CI runs on **GitHub Actions**. Two concerns: run every test surface
defined in [§5a](#5a-tests), and publish the bespoke container images
listed in [§2.2](#22-custom-logic-components) so the workflow has
something to `docker pull`.

**Workflows** (one YAML file each, under `.github/workflows/`):

- **`tests.yml`** — runs on every push and PR to any branch.
  1. `nf-core lint` (or `nextflow lint main.nf` if the nf-core scaffold
     is skipped — see [tasks/1_scaffold.md §9 q1](../tasks/1_scaffold.md)).
  2. **Container-coverage check** — assert every process in
     `modules/local/` resolves to a `container` directive via
     `conf/containers.config`, and that no tag is `latest`
     ([§1a](#1a-engineering-constraints)). A process without a
     container fails CI.
  3. **Stageable-file check** — grep-based lint asserts no process
     `script:` block contains a bare `${params.<file>}` interpolation
     that isn't wrapped in `file(...)` or on the scalar-param allow-list
     ([§1a](#1a-engineering-constraints)).
  4. **Python** — `unittest` on `bin/*.py` + `wf-report-boilerplate/`
     with `flake8`; each component's tests run inside its own
     container image ([§2.2](#22-custom-logic-components)) so the
     runtime environment matches production exactly.
  5. **Nextflow** — `nf-test` per-module + workflow-level; the
     workflow-level test is the P0 `-profile test -stub` end-to-end
     pass. Container runtime enabled so image-pull failures surface
     here.
  6. **Coverage report** — per-component coverage published as a job
     summary; <100% on non-mocked lines fails the build.

- **`images.yml`** — rebuilds bespoke container images. Triggered
  **only** on push to `main` (or on a release tag) touching one of:
  - `scripts/requirements.txt` (or per-image `containers/<name>/requirements.txt`)
  - `scripts/Dockerfile` (or per-image `containers/<name>/Dockerfile`)
  - the workflow file itself.

  **Not triggered by Python source changes** — images bundle the
  interpreter + pinned dependencies only; the source under
  `scripts/` / `bin/` is delivered to the container at runtime (see
  [§2.2](#22-custom-logic-components) container pattern). Editing a
  Python file rebuilds nothing, ships nothing, and iterates in
  seconds.

  For each of the five bespoke images from
  [§2.2](#22-custom-logic-components) — `wf5/samplesheet`,
  `wf5/coverage-gate`, `wf5/bin-target`, `wf5/barcode-validate`,
  `wf5/report`:
  1. Build with `docker/build-push-action`.
  2. Tag with the short git SHA of the commit that touched deps + (on
     release) the semver tag.
  3. Push to **DockerHub** — free for public repos, need to provide DockerHub
     API credential for image push.
  4. Emit a build-provenance attestation (`actions/attest-build-provenance`)
     — biosecurity reporting benefits from a verifiable chain from
     source to image.
  Path filters are per-image so a change to only `wf5/coverage-gate`'s
  requirements doesn't rebuild the other four.

- **`ref-build.yml`** *(optional, added when the reference-build script
  in [§4.4](#44-consolidated-build-script) is wired)* — scheduled
  monthly, invokes `scripts/build_refs.sh` in the `wf5/ref-build` image,
  publishes the resulting `refs/v<YYYY.MM>/` bundle to a release
  artefact (or object store). Manual `workflow_dispatch` trigger for
  ad-hoc rebuilds.

**Branch protection:** `main` requires `tests.yml` green before merge.
`images.yml` runs post-merge so the image tag matches the merged
commit.

**Local re-run:** every check in `tests.yml` must be runnable outside
CI (`nf-test`, `python -m unittest`, `flake8`, and the two grep
checks). No CI-only tooling.

## 6. Development phases

| Phase | Goal | Exit criteria |
|-------|------|---------------|
| **P0 — Scaffold** | Nextflow DSL2 skeleton, container plan, params schema, CI lint. | `nextflow run main.nf -profile test` completes a no-op end-to-end. |
| **P1 — Read prep + recruit + coverage gate** | Stages 1–6 wired with real test data. | Recruitment yield verified on plant test data at skim depth; coverage gate correctly subsamples an over-covered sample and soft-fails an under-covered one without blocking sibling samples. |
| **P2 — Assemble + polish + viz** | Stages 7–9. | Plant chloroplast assembled and visualised end-to-end. |
| **P3 — Bin + validate + extract** | Stages 10, 11, 13. | Barcode FASTA produced for plant + animal clean datasets; ORF validation passes. |
| **P4 — Annotate + collate + report** | Stages 12, 14, 15, 16; full output bundle, per-sample `report.html`, run-level `run-report.html`. | Taxodactyl handoff bundle validates against schema; full annotation produced; both report tiers render on multi-sample test invocation, including a mixed batch with one soft-failed sample. |
| **P5 — Contamination tolerance** | Minor-contamination datasets. | Acceptance criteria in §5 met. |
| **P6 — Negative-result clarity** | No-assembly and no-locus paths explicit in outputs. | Degraded-sample test yields explicit no-recovery diagnostic, distinct from null output ([brief.md §3.8](brief.md)). |

## 6a. Report design (per-sample `report.html`)

**Design reference:** [`mock-report.html`](mock-report.html). Use its **visual
layout, section structure, and interaction patterns** as the template.
**Do not** copy its scientific spec — the mock was produced for a SPAdes /
BUSCO-based pipeline, and its tool choices (assembler, completeness metric,
BLAST target DB, gene panel) must not leak into our specification. Content
comes from *our* pipeline stages (§2); the mock only informs how it's
presented.

### 6a.1 Design elements to adopt

- **Bootstrap-based single-file HTML** — self-contained, no external asset fetch, viewable offline.
- **Two-column header row:** *Key findings* (left, bulleted plain-English summary) + *Input parameters* (right, sample_id / kingdom / submitter metadata table).
- **Modal dialogs for detail views** — e.g. full parameter dump, tool versions with commit/build hashes. Keeps the main flow terse; detail one click away.
- **Multiple H2 top-level sections**, each a distinct part of the story, ordered as the pipeline runs.
- **Interactive charts** where a plot beats a table — read-length distribution, coverage-along-genome, per-contig coverage. Plotly is the mock's choice; keep unless it turns out to bloat the file substantially at multi-sample scale.
- **Inline SVG for the annotated organelle map** — clickable/hoverable gene features with tooltip metadata (name, strand, coordinates, source).
- **"Negative clarity" front-and-centre** ([brief.md §3.8](brief.md)) — when no organelle is assembled or no locus is extractable, the *Key findings* block says so explicitly at the top of the report, not buried in section 6.

### 6a.2 Section outline (our pipeline's content)

| H2 section | Content source | Notes |
|---|---|---|
| Header (H1 + two-column summary) | Sample metadata + top-line results | Key findings bullets: assembly outcome, coverage, top BLAST hit, gene-feature count, panel-locus extraction count. Right column: sample_id, kingdom, and any submitter-supplied optional columns from `samples.csv` (§0). |
| Sequencing data quality | `NANOPLOT_RAW` + `NANOPLOT_CLEAN` | Pre/post read stats, length + Q-score distributions, filter-yield summary from CHOPPER + FILTLONG. |
| Recruitment + coverage gate | `RECRUIT` + `COVERAGE_GATE` | Recruited-read count and total bases, estimated coverage vs configured MIN/MAX (§2.1), gate decision (pass / subsampled / soft-fail) with the seqtk fraction and seed if applicable. When soft-failed, this section is the terminal content — no assembly attempted. |
| *De novo* organelle assembly | `METAFLYE` + `MEDAKA` | Assembly statistics (contigs, N50, total bp), polishing summary, per-contig breakdown (length, coverage, circularised y/n). |
| Assembly quality assessment | `BANDAGE_NG` + `BIN_TARGET` | Assembly graph image, per-contig coverage-along-genome plot, binning classification (target / secondary / off-target). **No BUSCO** — mock includes it, our spec does not (kingdom-agnostic organelle BUSCO sets are patchy; revisit if needed). |
| Homology to reference databases | `BLAST_VALIDATE` | Top hits per target contig against the kingdom-appropriate RefSeq organelle DB (§4.2). Identity, coverage, subject accession. Explicit note when top-hit identity is below species-threshold — flags novel/underrepresented taxa without pretending to *assign* them (assignment is Taxodactyl's job — see [brief.md §2](brief.md) boundary). |
| Extracted barcode panel | `MINIPROT_EXTRACT` | Per-locus: extracted length, ORF status, genetic code used, coordinates, protein-identity to the miniprot reference. Panel loci not recovered are listed with the reason (not present / ORF broken / below length threshold). |
| Annotated organelle genome map | `ANNOTATE` | Inline SVG. Features from MITOS2 (animal) or the plant annotator TBD ([§8](#8-remaining-open-questions)); tooltip shows name, product, strand, coordinates, source score. |

### 6a.3 Elements to leave in the mock

- **SPAdes references / SPAdes-specific stats.** We use metaFlye — swap all assembler mentions.
- **BUSCO section.** Not in our spec. Do not add.
- **`core_nt` / their specific BLAST subject DB.** We validate against kingdom-partitioned RefSeq organelle DBs (§4.2).
- **Their proposed-taxonomy panel.** Downstream of our boundary — Taxodactyl produces taxonomy; our report reports homology hits only.
- **Their fixed gene panel details.** Loci come from Taxodactyl's accepted-loci config at runtime ([brief.md §3.7](brief.md)).

### 6a.4 Run-level `run-report.html`

Uses the same visual language as the per-sample report but content is a
cross-sample overview:

- Sample count + status breakdown: `ok` / `no_recovery` (assembled but no locus) / `low_coverage` (soft-failed at COVERAGE_GATE) / `hard_failure` (unexpected tool crash).
- Per-sample summary table — one row per sample_id, columns: kingdom, gate status, assembly outcome, coverage, top BLAST hit, panel recovery count. Each row links into `<sample_id>/report.html`.
- Run-level provenance panel: reference bundle version (§4.4 `manifest.json`), pipeline commit, samplesheet path + hash, invocation timestamp.

### 6a.5 Report implementation boilerplate

The [`wf-report-boilerplate/`](wf-report-boilerplate/) module is the
starting point for both report tiers. It is trimmed from a prior
Nextflow workflow to a clean scaffold that preserves the paradigms we
want to reuse:

- **Python + Jinja2 templating** with a single `render()` entry point
  ([`report.py`](wf-report-boilerplate/report.py)) invoked by the
  `COLLATE` / `RUN_REPORT` processes.
- **Self-contained single-file HTML** — all CSS, JS, and images inlined
  at render time via `_get_static_file_contents()`. No external asset
  fetch; report is viewable offline and archivable as-is.
- **Config-driven result-file discovery** ([`config.py`](wf-report-boilerplate/config.py)):
  glob-based `@property` accessors on `Config`, keyed off a per-sample
  `RESULT_DIR`. One property per output artefact; add as stages emit
  new files.
- **Class-based result objects** ([`results.py`](wf-report-boilerplate/results.py)):
  `AbstractDataRow` for single-row typed structs (e.g. `Metadata`,
  `RunQC`) and `AbstractResultRows` for tabular results driven by a
  CSV schema in `schema/`. Keeps template logic dumb; casting,
  bootstrap-class assignment, and derived fields live on the class.
- **Bootstrap 5 + Plotly + DataTables** shipped as vendored static
  assets; tabbed layout (`content-tabs.html`), stage-per-tab, with a
  reusable analyst-comment macro (`subjective.html`) and a
  save-as-read-only workflow (`save-modal.html` + `save-report.js`).

The boilerplate ships stubbed — concrete result classes, config
accessors, and stage-per-tab component templates are added as each
pipeline stage (§2) lands. [`mock-report.html`](wf-report-boilerplate/mock-report.html)
is retained inside the module as the visual design reference described
in §6a.

## 7. Reconciliation TODO against brief.md

Decisions made in this plan that supersede or refine brief.md v0.3:

- **§6 QC_FILTER → split into NANOPLOT_RAW + CHOPPER + FILTLONG + NANOPLOT_CLEAN.** Brief lists `chopper, NanoPlot`; this plan splits NanoPlot into pre/post passes and adds Filtlong for identity-weighted top-quality selection after chopper's threshold cut. Dorado handles adapter trimming at basecalling, so no separate trim stage is needed.
- **§6 add BANDAGE_NG, BLAST_VALIDATE, ORGANELLE_MAP.** Brief does not list visualisation or BLAST validation stages (ANNOTATE is now in brief.md §6).

## 8. Remaining open questions

Structural / architectural questions that need answering before or during
early development, distinct from the parameter tuning in §9. Bigger-picture
"which shape is the pipeline?" items land here; "what value does knob X
take?" items land in §9.

1. **Read-level fallback arm** for samples too degraded to assemble ([brief.md §8.2](brief.md)). Defer to post-P5 unless P5 testing surfaces a clear need. More likely to bite at skim depth.
2. **Assembly strategy — bait-then-refine?** Do we try a quick-and-dirty draft assembly with read baiting, and then re-map clean reads back to that to pick up any reads that baiting didn't capture? Then assemble those reads properly? Structural question (adds a whole pass), not a knob — parking here rather than in §9.
3. **Plant organelle annotator (stage 12 `ANNOTATE`, plant branch).** Tool choice deferred. Candidates: **[Chloe](https://github.com/ian-small/chloe)** (Julia, plastid-focused), **PGA**, **Mitofinder** (mt-only), or an in-house HMM-based annotator built on Prokka + custom HMM libraries. Note: any single tool may not cover both cp and mt; a split choice is on the table. Decision needed before P4. See [this paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9503105/pdf/ijms-23-10804.pdf) for deep dive.
4. **Organelle map renderer (stage 14 `ORGANELLE_MAP`).** Tool choice deferred. Candidates: **[pyCirclize](https://github.com/moshi4/pyCirclize)** (Python, actively maintained, GenBank input), **[plotgenes](https://github.com/aaronphillips7493/plotgenes)** (from CLAW's author), Circos with a custom GenBank→config converter, or an in-house Jinja+SVG renderer that reads GFF/GenBank directly (would fit the report renderer container per [§6a.5](#6a5-report-implementation-boilerplate)). Inline SVG output preferred so it drops into the report. Decision needed before P4.

## 9. Fine tuning (post-prototype benchmarking)

Items to benchmark **after** the first end-to-end prototype is running on
real data (roughly end of P3). Each entry is a knob or design choice where
we've picked a sensible default from prior art or first principles, but the
right value / choice will only become obvious with our own test data in
hand.

**Format:** `knob` → `current default` → `experiment` → `success metric`.

| # | Knob / choice | Current default | Experiment | Success metric | Cross-ref |
|---|---|---|---|---|---|
| 1 | **RECRUIT filter strictness** (plant cp) | `samtools view -F 4 -q 1` (any primary mapping) | Compare against ptGAUL-style PAF filter: `alignment_length/read_length ≥ 0.7` and `alignment_length ≥ 1 kb`. | Cleaner input to METAFLYE (fewer off-target reads) without a recall hit that pushes COVERAGE_GATE below MIN on borderline samples. Measure: recruited-read count, estimated coverage, final assembly N50 + edge count. | §2 stage 5, §3.6 |
| 2 | **COVERAGE_GATE MIN/MAX** (per kingdom / organelle) | plant cp 30×/500×; animal mt 30×/300× | Sweep MIN downward on marginal plant / animal samples to find the floor where assembly still recovers a canonical structure. Sweep MAX to find the point where metaFlye graph tangles reappear. | Widest MIN–MAX band that still produces canonical assemblies. Failed samples at MIN vs. successful ones just above tells us the real floor. | §2.1.1 |
| 3 | **Flye `--genome-size` hint** | plant `2m`, animal `20k` | Vary the plant hint (500 kb, 1 m, 2 m, 4 m) on a plant sample where mt is present. Does under-hinting fragment mt? Does over-hinting hurt cp assembly? | Assembly contiguity for whichever organelle is under scrutiny (N50, edge count). | §3.4 |
| 4 | **Flye `--asm-coverage`** (plant) | `MAX − 20` (= 480 with default MAX 500) | Vary in 50× steps around 480; also try without the flag. | Plastid edge-count + total assembly length converging on canonical 3-edge quadripartite. | §3.5 |
| 5 | **Flye read-mode** (`--nano-hq` vs `--nano-corr` after upstream correction) | `--nano-hq` on raw Dorado SUP R10.4.1 | Insert an error-correction stage (HERRO / Canu-correct / Ratatosk) between FILTLONG and RECRUIT; run Flye with `--nano-corr`. Compare against the default. | Assembly indel rate + panel ORF pass rate. If gain is modest (<5% ORF-pass improvement) the runtime cost isn't justified. | §2 stage 7, §8 (was open Q4) |
| 6 | **Polish depth** (medaka iterations) | `--polish` opt-in, single pass when enabled | For samples that fail ORF validation without polish, compare: no polish vs. single medaka pass vs. iterative (2–3 passes). | Panel ORF pass rate at each polish depth vs. wall-clock cost. Look for diminishing returns. | §2 stage 8, brief.md §8.3 |
| 7 | **CHOPPER thresholds** | min length + min mean Q from config (values TBD in P1) | Sweep min length (500 / 1000 / 2000 bp) and min Q (7 / 10 / 12) on plant and animal skim data. | Recruited-read yield after RECRUIT and estimated coverage post-gate. Aim for the loosest filter that doesn't push borderline samples below COVERAGE_GATE MIN. | §2 stage 2 |
| 8 | **FILTLONG `--keep_percent`** | 95 | Compare 90 / 95 / 99 / disabled. FILTLONG removes reads by identity-weighted quality percentile — how much low-Q tail does metaFlye actually tolerate at skim depth? | Assembly contiguity + panel ORF pass rate vs. recruited-read count. | §2 stage 3 |
| 9 | **Assembler choice** ([brief.md §8.1](brief.md)) | metaFlye | Benchmark against GetOrganelle, Oatk, and (for plants) ptGAUL itself as an end-to-end alternative. | Assembly completeness (canonical structure, panel gene recovery), indel rate after polish, within-kingdom mixture behaviour, runtime, packaging complexity. | §2 stage 7, brief.md §8.1 |
| 10 | **BIN_TARGET coverage-spike threshold** | `2× the median non-target coverage` (working assumption) | Once real data exists, plot per-contig coverage distributions and pick a threshold empirically. Currently a first-principles guess. | Correct target-contig selection rate on the P3 clean test datasets (all four kingdom-organelle combinations). | §2 stage 10 |
| 11 | **Miniprot per-locus representative count** | 5–10 proteins per locus per kingdom | Vary between 1 (single canonical protein), 5–10 (current), and 50+ (broad clade coverage). Measure sensitivity on divergent test taxa. | Panel locus recovery rate on the most divergent test sample per kingdom. Runtime scaling for the top end. | §4.3 |

Each entry becomes a P-something task once the prototype is stable enough
to run these experiments repeatably. Grouping suggestion:

- **P-tune-A (read prep + recruit):** items 1, 7, 8. All touch stages 2–5 and can share test data.
- **P-tune-B (assembly parameters):** items 2, 3, 4, 5, 6, 9. Shared assembly benchmark harness.
- **P-tune-C (downstream):** items 10, 11. Each is stage-local.

Priority ranking will fall out of what breaks in P3–P5 testing; nothing
here blocks P0–P4.
