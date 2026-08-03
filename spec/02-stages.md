## 2. Stage detail

| # | Stage | Tool | Inputs | Outputs | Ref data | Key params / notes |
|---|-------|------|--------|---------|----------|--------------------|
| 1 | `NANOPLOT_RAW` | NanoPlot | raw FASTQ | read length/quality plots, summary TSV | — | Baseline metrics; feeds the report. Adapter/barcode trimming is performed upstream by Dorado at basecalling time, so no separate trim stage is needed. |
| 2 | `CHOPPER` | chopper | raw FASTQ | length-/quality-filtered FASTQ | — | Min length and min mean Q from config. Minimal filtering (preserve target reads). ONT-native. |
| 3 | `FILTLONG` | Filtlong | chopper-filtered FASTQ | identity-weighted top-quality FASTQ | — | Complements chopper's threshold cut with a percentile-based, identity-weighted selection (Filtlong scores reads by window quality and length). Prevents low-identity long reads from swamping METAFLYE at skim depth. Config: `--min_length`, `--keep_percent` (default 95). |
| 4 | `NANOPLOT_CLEAN` | NanoPlot | filtlong output | post-clean read stats | — | Pre/post comparison goes into report. |
| 5 | `RECRUIT` | minimap2 + samtools + seqtk | cleaned FASTQ, single organelle reference (`${assembly_target}.mmi`) | `recruited.fastq.gz` | [§4.1](04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle) | Positive recruitment ([brief.md §3.2](brief.md)). One index per run, selected by `meta.assembly_target` per [§1a](01-pipeline-flow.md#1a-engineering-constraints). Command shape: `minimap2 -ax map-ont -t $T $refs/${assembly_target}.mmi $reads \| samtools view -b -F 4 -q 1 -@ $T` — extract mapped read IDs, then `seqtk subseq` back to FASTQ. Liberal threshold — coarse enrichment only. Non-recruited reads discarded ([brief.md §3.3](brief.md)). Pattern adapted from [CLAW](../CLAW/Snakefile). **Optional stricter filter** (P1 benchmark, `plant_pt` candidate): PAF-based recruitment with `minimap2 -cx map-ont`, retaining reads where alignment length / read length ≥ 0.7 and alignment length ≥ 1 kb — the [ptGAUL](../reference-material/ptgaul/ptGAUL.sh) filter. Trades a small recall hit for cleaner input to METAFLYE; worth evaluating on `plant_pt` samples where plastid recruitment is high-signal. |
| 6 | `COVERAGE_GATE` | seqkit stats + minimap2 + seqtk sample | recruited FASTQ, organelle reference panels, per-target coverage limits (§2.1) | `gated.fastq.gz` (subsampled if needed), `coverage.json`, sample-status marker | [§4.1](04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle) | Sum recruited bases via `seqkit stats -T`; where the target has a sibling organelle panel, assign each read to its best-matching panel and estimate coverage as `target_assigned_bases / nominal_organelle_size` (size fixed per `assembly_target`) — see [§2.1.5](#215-sibling-organelle-carry-over-in-the-estimate). If `< MIN` → **soft-fail** the sample: emit no-recovery marker, skip downstream, sample proceeds to `COLLATE` with an explicit low-coverage `report.html` (§2.1.3). If `> MAX` → `seqtk sample -s $seed $recruited $frac` where `frac = MAX / estimated`. If in-band → passthrough. Cross-sample isolation: `errorStrategy 'ignore'` on this process; failure state is data (a marker file), not a Nextflow error, so other samples are unaffected. |
| 7 | `METAFLYE` | Flye `--meta` (`--nano-hq` \| `--nano-corr`) | gated FASTQ | `assembly.fasta`, `assembly_graph.gfa`, `assembly_info.txt` | — | Retain `assembly_info.txt` — per-contig coverage is load-bearing for `BIN_TARGET`. `--meta` retained to tolerate low-level contamination ([brief.md §3.5](brief.md)). **Read-mode flag depends on whether an error-correction stage precedes Flye** — see §9 item 5. Default assumption: raw Dorado SUP → `--nano-hq`. Pass an organelle-keyed `--genome-size` hint (see §3.4). Note: `--asm-coverage` is incompatible with `--meta` in Flye 2.9.6 and is not used. **Polishing:** rely on Flye's built-in polisher (`--iterations 1`, the Flye default for `--nano-hq`). A separate `MEDAKA` stage is deferred — see §9 item 6. See the [Flye user guide](../reference-material/flye-user-guide.md) for parameter reference when tuning any of the Flye knobs in §9. |
| 8 | ~~`MEDAKA`~~ | — | — | — | — | **Deferred.** Original plan was an opt-in `medaka_consensus` polish after METAFLYE. Shelved pending benchmark evidence that Flye's built-in iteration is insufficient for ORF recovery — tracked as §9 item 6. Stage 8 kept as a numbered slot so downstream stage references remain stable. |
| 9 | `BANDAGE_NG` | BandageNG | `assembly_graph.gfa` | graph image (PNG/SVG) | — | Diagnostic only; helps human review of tangled assemblies / multi-organism cases. |
| 10 | `BIN_TARGET` | minimap2 + custom | assembly FASTA, `assembly_graph.gfa`, `assembly_info.txt`, declared organelle ref (`${assembly_target}.mmi`) **+ sibling organelle ref(s)** | target contig(s) FASTA, secondaries TSV, plastid isoforms (`plant_pt` only) | [§4.1](04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle) | Per-contig binning per [§3.7](03-organelles.md#37-target-binning-criteria-per-assembly-target): merged-interval aligned fraction ∩ identity against the declared panel ∩ declared panel beating every sibling panel. **cp/mt discrimination inside the sample is required** on plant targets — plastid carry-over dominates `plant_mt` assemblies by coverage ([§3.3](03-organelles.md#33-specific-issues-and-decisions), [§3.7.1](03-organelles.md#371-homology-is-measured-against-sibling-panels-too)). Coverage ranks admitted candidates for dominance; it is **not** an admission gate ([§3.7.3](03-organelles.md#373-coverage-is-a-ranking-signal-not-a-gate)). Secondaries logged not emitted ([brief.md §3.5](brief.md)). **`plant_pt` branch:** apply the quadripartite canonicalisation algorithm from §3.6 to emit `path1.fasta` (primary) + `path2.fasta` (alternative isoform), subject to the §3.6 step 5 selection precondition; record edge-count structural signal in report. |
| 11 | `BLAST_VALIDATE` | blastn | target contigs, target-specific BLAST DB | per-contig top hits + identity/coverage TSV | [§4.2](04-reference-data.md#42-refseq-organelle-dbs-for-blast_validate) | Sanity check that selected contigs match the declared `assembly_target`. BLAST DB fixed per-run by `meta.assembly_target`. Feeds metadata + report. |
| 12 | `ANNOTATE` | MITOS2 (animal mt) / **TBD** (plant cp+mt — see [§8](07-open-questions.md#8-remaining-open-questions)) | target contigs | annotation GFF + GenBank | — | Full organelle annotation ([brief.md §3.6](brief.md)). Feeds `ORGANELLE_MAP` and provides annotated organelle as a standalone output. Kingdom-dispatched: animal → MITOS2; plant → tool TBD. |
| 13 | `MINIPROT_EXTRACT` | miniprot | target contigs, locus panel (protein refs subset for the sample's `assembly_target`) | barcode FASTA, coordinates GFF, per-locus validation | [§4.3](04-reference-data.md#43-protein-panel-for-miniprot_extract) | Protein-to-genome; works on divergent taxa with no same-species reference. Validates length, ORF, target-appropriate genetic code, internal stops ([brief.md §7](brief.md)). Locus panel is subset per `meta.assembly_target` (rbcL/matK for `plant_pt`; cox1/nad1 for `plant_mt`; COX1/CYTB for `animal_mt`). Runs independently of ANNOTATE — the barcode panel is the contractual output; ANNOTATE is supplementary. |
| 14 | `ORGANELLE_MAP` | **TBD** (see [§8](07-open-questions.md#8-remaining-open-questions)) | annotated GenBank from `ANNOTATE` | organelle map (PNG/SVG) | — | Diagnostic visualisation from the full annotation. Renders a circular/linear feature map for both plant (cp, mt) and animal (mt) targets. Tool choice deferred. |
| 15 | `COLLATE` | — | all upstream outputs (per sample) | `outdir/<sample_id>/{organelle_assembly.fasta, organelle_annotation.gff, barcodes.fasta, metadata.json, report.html}` | — | Per-sample aggregation; produces the Taxodactyl handoff bundle ([brief.md §5](brief.md)). Handles soft-failed samples: if the `COVERAGE_GATE` marker says fail, emit a minimal bundle with a low-coverage `report.html` and no assembly outputs. |
| 16 | `RUN_REPORT` | — | all `COLLATE` outputs + samplesheet | `outdir/run_manifest.json`, `outdir/run-report.html` | — | Cross-sample summary and manifest. Records samplesheet snapshot, reference bundle version (§4.4 `manifest.json`), pipeline commit, per-sample success/no-recovery/low-coverage breakdown. |

## 2.1 Coverage gate (§2 stage 6)

Post-recruitment, before assembly, estimate coverage from recruited bases and
either subsample or soft-fail. Purpose: (a) keep metaFlye out of the
over-coverage regime where the graph tangles, and (b) fail fast when there is
too little signal to recover an organelle, without letting a bad sample block
a batch.

### 2.1.1 Per-target limits

Config-driven table (default values below; overridable via
`params.coverage_limits.<assembly_target>`):

| Assembly target | Nominal size (for coverage denominator) | MIN cov | MAX cov | Notes |
|---|---|---|---|---|
| `plant_pt` | 150 kb | 30× | 500× | Plastid — small, very high copy. MAX generous but capped to keep the metaFlye graph tractable. |
| `plant_mt` | 400 kb | 30× | 300× | Plant mitogenome — larger, medium copy, recombinant structure. Coverage more variable across taxa; MAX slightly tighter. |
| `animal_mt` | 17 kb | 30× | 300× | Tight organelle, high per-cell copy — MAX is generous because Flye handles small circles well up to a few hundred ×. |

One row per target flows through the pipeline independently, so limits
apply directly per run without any cross-organelle gating logic.

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

### 2.1.5 Sibling-organelle carry-over in the estimate

**The limitation.** [§2.1.2](#212-estimation-formula) totals *every*
recruited base against a nominal size fixed per `assembly_target`. RECRUIT
is coarse positive enrichment and deliberately does not separate the two
plant organelles ([CONSTITUTION principle 5](../CONSTITUTION.md)), so on
`plant_mt` the recruited pool carries a large plastid fraction — the
abundance gradient runs plastid → mitochondrion in plant tissue, the two
share real homology, and the `plant_mt` panel itself contains
plastid-derived (NUPT) sequence. Estimating from all recruited bases then
over-states mitochondrial depth, and the gate passes samples that should
soft-fail. This is [principle 7](../CONSTITUTION.md) in reverse: the gate
exists to stop a degraded sample producing confident-looking output, and
its own input was what hid the degradation.

Measured on the `INT-PLANT-01-mt` fixture (task 25 §2): **70.5 %** of
recruited bases were plastid, the gate reported 92.60× and `status: "ok"`,
and the mitochondrial-only estimate was **27.29×** — below the 30× MIN.
The assembly that followed was 61 % chloroplast by length, with the
genuine mitochondrial contigs at 7–8×.

**The resolution.** Where the declared target has a sibling organelle
panel, C2 maps the recruited pool against the declared panel and every
sibling, assigns each read to the panel it aligns to best on the
merged-interval metric of
[§3.7.2](03-organelles.md#372-aligned-fraction-is-merged-across-all-alignment-blocks),
and gates on target-assigned bases only:

```
estimated_cov = target_assigned_bases / nominal_organelle_size
```

Sibling map is the one C3 applies
([§3.7.1](03-organelles.md#371-homology-is-measured-against-sibling-panels-too)):
`plant_mt` ↔ `plant_pt`; `animal_mt` has no sibling and is unaffected.
Ties fall to the sibling — a read equally explained by both organelles is
not evidence of target depth ([rule 18](../CONSTITUTION.md)). Subsampling
is unchanged: `seqtk sample` still draws uniformly from the whole pool, so
the realised reduction applies to the target-assigned fraction too.

`sample_status.json` records `coverage_basis`, `total_recruited_bases`,
`target_assigned_bases`, `sibling_assigned_bases`,
`sibling_organelle_fraction`, and `sibling_panels_scored`, so an auditor
can see both the corrected estimate and what it corrected for.

**Fallback.** Where a required `.mmi` is absent from the bundle or
minimap2 fails, C2 warns to stderr, records `coverage_basis:
"total_recruited"` with `sibling_panels_scored: []`, and estimates from
the whole pool as before. A degraded reference bundle must not fail the
sample ([principle 8](../CONSTITUTION.md)).

**Post-assembly cross-check.** C3 additionally reports
`sibling_carryover` in `bin_metadata.json` — the fraction of *assembled*
bases binned as `sibling_organelle`, with a `sibling_carryover_warning`
above `params.bin_target_thresholds.<target>.sibling_warn_fraction`
(default 0.30, provisional pending
[§9 item 10](07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)).
This is the same signal measured where it is unambiguous, and it surfaces
in the report ([§6a.2](06a-reports.md#6a2-section-outline-our-pipelines-content)).

**Not resolved here.** Filtering the sibling reads out at RECRUIT would
fix the root cause and stop METAFLYE spending its depth budget on the
wrong organelle. That sits on a [principle 5](../CONSTITUTION.md)
boundary and risks discarding genuine NUPT-spanning mitochondrial reads,
so it is held for
[§9 item 1](07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
(RECRUIT filter strictness) rather than pre-empted.

## 2.2 Custom logic components

Most stages are thin wrappers around off-the-shelf bioinformatics tools
(`minimap2`, `Flye`, `blastn`, `miniprot`, etc.) invoked
directly from a `script:` block. This section catalogues the places
where custom logic — Python scripts, non-trivial parsing, or novel
decision code — is required. Each item is:

- an in-house script living under `bin/` (Nextflow auto-stages `bin/`
  onto every process `PATH`);
- distributed inside a **small dedicated container** (`containers/<name>/`)
  per [§1a](01-pipeline-flow.md#1a-engineering-constraints) — no host-Python fallback;
- covered by unit tests independent of the Nextflow harness, so the
  logic is testable without spinning up the pipeline.

If you find yourself writing more than ~10 lines of `awk`/`bash` inside
a `.nf` script, promote it into `bin/` and add an entry here.

| # | Component | Container | Stage(s) | Purpose | Spec anchor |
|---|---|---|---|---|---|
| C1 | `parse_samplesheet.py` | `wf5/samplesheet:<tag>` (thin Python + `pandas`/`csv` + `pyyaml`) | [§2 stage 0](#2-stage-detail) | CSV validation beyond nf-schema's reach: `sample_id` uniqueness across the sheet, pipe-split of `reads`, absolute-path rejection with an actionable error, resolution against `--data-dir`, existence check of every resolved path, `assembly_target` normalisation + enum check `{animal_mt, plant_pt, plant_mt}`, ISO 8601 date warning-not-fail. Emits a normalised JSON array that Nextflow `splitJson`-consumes into the per-sample channel. | [§0](00-overview.md#0-input-sample-sheet) |
| C2 | `coverage_gate.py` | `wf5/coverage-gate:<tag>` (Python stdlib only + `seqkit`/`seqtk`/`minimap2` in `PATH`) | [§2 stage 6](#2-stage-detail) | Reads `seqkit stats -T` output, splits the recruited pool by organelle panel where the target has a sibling ([§2.1.5](#215-sibling-organelle-carry-over-in-the-estimate)), computes `estimated_cov = target_assigned_bases / nominal_organelle_size` from the per-target limits table, executes the three-branch decision (soft-fail / passthrough / subsample), invokes `seqtk sample -s $seed` when subsampling, and writes `sample_status.json` + `coverage.json`. **Always exits 0** so a coverage decision is data, not a Nextflow error ([§2.1.4](#214-cross-sample-failure-isolation)); a missing sibling index degrades to a whole-pool estimate rather than failing the sample. `minimap2` is in the image from task 25 — C2 is no longer stdlib+`seqkit`/`seqtk` only, and now takes the reference bundle as a staged input. | [§2.1](#21-coverage-gate-2-stage-6) |
| C3 | `bin_target.py` | `wf5/bin-target:<tag>` (Python + `biopython` + `pysam`/`mappy` for GFA/BAM parsing) | [§2 stage 10](#2-stage-detail) | Per-contig classification per [§3.7](03-organelles.md#37-target-binning-criteria-per-assembly-target): (a) merged-interval aligned fraction + identity vs the declared panel (`${assembly_target}.mmi`), (b) the declared panel out-scoring every sibling organelle panel, with (c) coverage as the dominance ranking among admitted candidates and (d) longest-ORF as a recorded diagnostic only. Thresholds are per-target ([§3.7.6](03-organelles.md#376-revised-per-target-criteria)). Selects the dominant target contig (all candidates for `plant_mt`), records secondaries in a diagnostic TSV, and records circularity from Flye's `circ.` column with end-overlap fallback ([§3.7.4](03-organelles.md#374-circularity-comes-from-flye-not-from-end-overlap)). | [§3.3](03-organelles.md#33-specific-issues-and-decisions), [§3.7](03-organelles.md#37-target-binning-criteria-per-assembly-target), [§9 item 10](07-open-questions.md#9-fine-tuning-post-prototype-benchmarking) |
| C4 | `plastid_canonicalise.py` — in-house implementation of the algorithm in [§3.6](03-organelles.md#36-plastid-quadripartite-canonicalisation), specified in full by [plastid-canonicalisation.md](plastid-canonicalisation.md). | shares `wf5/bin-target:<tag>` (invoked by C3 on the plant-cp branch) | [§2 stage 10](#2-stage-detail), plant-cp branch only | Parses Flye's `assembly_graph.gfa`, classifies edges by length + depth (LSC = longest, IR = deepest, SSC = remainder for the 3-edge case), emits `path1.fasta` and `path2.fasta` for the two SSC-orientation isoforms. Handles the 1-edge (resolved circle passthrough) and N≠3 (diagnostic-only) cases explicitly. Written in-house from the algorithm specification. | [§3.6](03-organelles.md#36-plastid-quadripartite-canonicalisation) |
| C5 | `validate_barcodes.py` | `wf5/barcode-validate:<tag>` (Python + `biopython`) | [§2 stage 13](#2-stage-detail) | Post-processes miniprot output: per-locus length check, ORF check under the target-appropriate NCBI genetic-code table (table 11 for `plant_pt`, table 1 for `plant_mt`, table 2/5 clade-trial for `animal_mt`), internal-stop check, protein-identity threshold. Emits `barcodes.fasta` (validated only) + a per-locus `validation.tsv` recording pass/fail reasons for panel loci that dropped out. **`animal_mt` clade trial** ([§3.3](03-organelles.md#33-specific-issues-and-decisions)): attempts NCBI tables 2 (vertebrate) and 5 (invertebrate), picks the one with a valid ORF, records the chosen table in metadata. | [§3.2](03-organelles.md#32-per-stage-parameter-selection), [§3.3](03-organelles.md#33-specific-issues-and-decisions) |
| C6 | `collate.py` | `wf5/report:<tag>` (shared with C7 and the report renderer — Python + Jinja2 per [§6a.5](06a-reports.md#6a5-report-implementation-boilerplate)) | [§2 stage 15](#2-stage-detail) | Detects the [`COVERAGE_GATE`](#21-coverage-gate-2-stage-6) soft-fail marker and dispatches to either the full or minimal per-sample bundle, aggregates upstream outputs into the layout in [§0](00-overview.md#0-input-sample-sheet), emits `metadata.json` (sample meta + tool versions + reference-bundle version), and invokes the [report renderer](../wf-report-boilerplate/report.py) to write `report.html`. | [§6a.5](06a-reports.md#6a5-report-implementation-boilerplate) |
| C7 | `run_report.py` | shares `wf5/report:<tag>` with C6 | [§2 stage 16](#2-stage-detail) | Cross-sample join: reads every per-sample `metadata.json`, classifies each sample into `ok` / `no_recovery` / `low_coverage` / `hard_failure`, emits `run_manifest.json` (samplesheet snapshot + reference-bundle version from [§4.4 manifest](04-reference-data.md#44-consolidated-build-script) + pipeline commit + invocation timestamp), and renders the run-level HTML via the same Jinja machinery as C6. | [§6a.4](06a-reports.md#6a4-run-level-run-reporthtml) |

**Non-workflow custom logic (pipeline-adjacent):**

- `scripts/build_refs.sh` (or `.py`) — the consolidated reference-bundle
  builder described in [§4.4](04-reference-data.md#44-consolidated-build-script). Not
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
off-the-shelf stages (see [§1a](01-pipeline-flow.md#1a-engineering-constraints)).

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
`Dockerfile` changes ([§5b](05-test-data.md#5b-ci)); day-to-day Python edits ship
via the next `nextflow run` without a Docker step.

**Test surface:** each `bin/*.py` gets a unit-test file under
`scripts/tests/` covering its decision logic against synthetic
fixtures, with external tool calls (`seqkit`, `seqtk`, `minimap2`,
etc.) mocked at the subprocess boundary so the module is imported
in-process rather than spawned — this is also what makes the module's
branches visible to `coverage` (see
[§5a](05-test-data.md#5a-tests)/[§5b](05-test-data.md#5b-ci)).
Nextflow-level integration testing is separate (P0 `-profile test`
covers channel wiring; P1+ tests exercise real tools). The unit tests
run in the shared `neoformit/daff-wf5-scripts:test` image in CI, not
per-component containers — see [§5b item 4](05-test-data.md#5b-ci).
