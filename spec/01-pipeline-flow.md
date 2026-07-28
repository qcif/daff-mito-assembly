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

- **One sample row → one organelle assembly.** Each row in the sample
  sheet declares an `assembly_target` from a fixed set:
  `animal_mt`, `plant_pt`, `plant_mt`. Every stage — recruitment,
  coverage gate, assembly, polish, binning, validation, annotation,
  extraction — operates against that single target. A plant sample
  whose operator wants both plastid and mitogenome assembled adds
  **two rows** (same reads, `assembly_target=plant_pt` on one,
  `assembly_target=plant_mt` on the other) and the two runs stay
  independent through the entire pipeline. This keeps reference
  selection, size hints, coverage limits, genetic-code tables, and
  BLAST DBs unambiguous per-run without any within-sample fan-out
  logic. See [spec §3.2](03-organelles.md#32-fork-vs-dynamic-parameter--by-stage)
  for the per-stage effect and [spec §0](00-overview.md#0-input-sample-sheet)
  for the samplesheet schema.



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
  [§4.4](04-reference-data.md#44-consolidated-build-script) may use either pattern. If in
  doubt, prefer pattern 1 — it's harder to misuse.
