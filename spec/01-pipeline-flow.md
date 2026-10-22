## 1. Pipeline flow

FastP - long - alternative for QC? Ont suite alternatives?

```
samples.csv
   │
   ▼
[0] PARSE_SAMPLESHEET ─────────► one channel item per sample: (meta, reads, intake.json)
   │                              meta = {sample_id, kingdom, sample_info}
   │                              validates gzip + FASTQ integrity of every input;
   │                              invalid reads ─► straight to COLLATE (`invalid_input`)
   │                              repeated read IDs: identical copies ─► DEDUPLICATE (below),
   │                              conflicting copies ─► `invalid_input`
   ▼    per-sample fan-out from here
[0b] DEDUPLICATE (conditional) ► one copy per read ID (seqkit rmdup); only when
   │                              intake.json sets dedup_required
   ▼
[1] NANOPLOT_RAW ──────────────► raw read stats (HTML/PNG), on the *submitted*
   │                              reads, duplicates included
   │                              (adapters already trimmed by Dorado)
   ▼
[2] CHOPPER ──────────────────► quality-/length-filtered reads (deduplicated
   │                              reads where DEDUPLICATE ran)
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
   │                              assemble-but-flag if < WARN;
   │                              fail sample (soft) if < HARD_MIN
   ▼
[7] METAFLYE (--meta --nano-hq)  draft organelle assembly + assembly graph
   │                              (Flye built-in polish, --iterations 1)
   ▼
[8] — reserved (MEDAKA deferred; see §7 open Q6)
   │
   ├──────────────────────────────┐
   ▼                              ▼
[9] BIN_TARGET (coverage + ref)  [10] BANDAGE_NG (+ C12) ──► assembly graph
   │   target contig(s) +           │   SVG, node-labelled, coloured by [9]'s
   │   secondaries                  │   binning bucket (diagnostic;
   │                                │   sentinel-colour round-trip recovers
   │                                │   the segment names Bandage's SVG
   │                                │   omits — joined to [9]'s
   │                                │   bin_metadata.json only at the final
   │                                │   annotate step, off the critical path)
   ▼
[9a] ROTATE_ANCHORS → ROTATE_TARGET (C13) ► circular `animal_mt` contigs
   │   reverse-complemented / rotated so COX1 starts at position 1;
   │   a circle that holds the genome twice is collapsed to one copy
   │   (other targets pass through; §3.8)
   ▼
[11] BLAST_VALIDATE ──────────► per-contig identity to kingdom organelle refs
   │
   ▼
[12] MINIPROT_CDS ────────────► all protein-coding features (cds.gff)
   │                              one broad miniprot pass over the
   │                              comprehensive protein panel (§4.3);
   │                              feeds BOTH [13] and [13a]
   ├──────────────┐
   ▼              ▼
[13] EXTRACT_   [13a] ANNOTATE ► organelle annotation (GFF3)
     BARCODES ──► barcode FASTA        + annotation_summary.json
   │              coords + ORF         CDS from [12]; tRNA/rRNA from
   │              validation.          MITOS2 (animal_mt) or barrnap
   │              Subset of [12]       (plant_pt, BARRNAP step first).
   │              by assets/loci.json  plant_mt is CDS-only (§8 Q3).
   │              by assets/loci.json
   ▼              ▼
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
  coverage gate, assembly, binning, validation, annotation,
  extraction — operates against that single target. A plant sample
  whose operator wants both plastid and mitogenome assembled adds
  **two rows** (same reads, `assembly_target=plant_pt` on one,
  `assembly_target=plant_mt` on the other) and the two runs stay
  independent through the entire pipeline. This keeps reference
  selection, size hints, coverage limits, genetic-code tables, and
  BLAST DBs unambiguous per-run without any within-sample fan-out
  logic. See [spec §3.2](03-organelles.md#32-per-stage-parameter-selection)
  for the per-stage effect and [spec §0](00-overview.md#0-input-sample-sheet)
  for the samplesheet schema.



- **Every task image provides a procps `ps` (not BusyBox).** The
  execution report, timeline and trace are always on, and Nextflow's
  metrics collection runs `ps` before the task script. A task whose image
  lacks it fails outright, and under `errorStrategy = 'ignore'` the
  sample then drops out silently.

- **A failing tool must fail its process.** Process scripts run under
  `bash -euo pipefail` (`process.shell` in `conf/base.config`), so a
  non-zero exit anywhere in a pipe fails the task. A pipe stage whose
  non-zero exit is expected must test that specific status, never
  `|| true`. `bin/` scripts that run subprocess pipelines must check every
  stage's return code (use `bash -o pipefail -c`, or `check=True` on each
  stage).

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
- **Prefer biocontainers or nf-core-published images.** Container
  selection order:
  1. **Single-tool biocontainer** on quay.io/biocontainers — the
     default for any stage wrapping one off-the-shelf tool.
  2. **`mulled-build` bundle** when a process needs multiple
     bioconda-packaged tools together (e.g. `RECRUIT` needs
     `minimap2` + `seqtk` + `samtools`). Do **not** hand-write a
     Dockerfile for this case — use `galaxy-tool-util`'s
     `mulled-build` to compose the bioconda deps into a
     content-hashed image and mirror it to our registry:
     ```sh
     pipx install galaxy-tool-util
     mulled-build build 'minimap2=2.31,seqtk=1.5,samtools=1.24'
     docker tag \
         quay.io/biocontainers/mulled-v2-<hash>:<build-hash> \
         neoformit/daff-wf5-<name>:<pinned-tag>
     docker push neoformit/daff-wf5-<name>:<pinned-tag>
     ```
     The hashed quay.io tag is deterministic given the same deps, so
     the bundle is reproducible without a Dockerfile in the repo.
     Record the exact `mulled-build` invocation in
     `tasks/2_containers.md` for provenance.
  3. **Hand-written Dockerfile** in `containers/<tool>/` only when no
     bioconda package exists or the image needs non-bioconda system
     tooling. Rare — reach for it last.
  Bespoke images for **custom Python logic** are separate
  ([spec §2.2](02-stages.md#22-custom-logic-components)) and follow
  their own "deps in the image, code at runtime" pattern.
- **Directory outputs are never empty.** Azure blob storage has no real
  directories: a "directory" is only a shared prefix of the files in it,
  so an empty one does not exist once uploaded, and a process declaring
  it as an output fails with "missing output" on Azure Batch while
  passing locally. Every `path("dir")` output that can be empty (in the
  real script or the stub) gets a `touch dir/.keep` in both `script:` and
  `stub:`. Use `.keep` everywhere. It is hidden, so `*.ext` globs skip it,
  and consumers must never match or count it. The lint job's
  directory-output check fails any module with a directory output that
  neither writes `dir/.keep` nor sits on its explicit allow-list (with
  a reason). Directories in the reference bundle follow the same rule.
- **File inputs must be stageable, not raw path strings.** On a remote
  executor (AWS Batch, Azure, Kubernetes, HPC scratch) Nextflow only
  copies a file to the compute node if it can see a Nextflow `Path`
  object referring to it. The way to produce one is:
  1. **Declare as an `input: path` channel.** Materialise
     the param into a value channel at the top of `main.nf`
     (`ch_blastdb = Channel.value(file(params.blastdb))`) and pass it as
     `input: path blastdb` on every consuming process. This makes the
     dependency explicit in the process signature, keeps `script:`
     blocks free of `params.*` references, and lets a process be tested
     in isolation with a fixture path.
  Do **not** rely on inline `${file(params.x)}` in a `script:` block:
  interpolation turns the `Path` into a string the executor cannot
  stage. It only works where the path is visible from the compute node
  (shared filesystem), and fails on Azure Batch (task 57: a `blast_db`
  referenced this way resolved to `/refdata-wf5/...` on the node).

  What **does not** work: bare `${params.blastdb}` in a script block —
  that interpolates to a plain string, and the executor has nothing to
  stage. Also never write to paths outside the process work directory,
  and don't read from `$projectDir`/`$baseDir` at runtime for data.

  Rule of thumb: per-sample data and any file that fans out always go
  through explicit channels; workflow-wide reference bundles listed in
  [§4.4](04-reference-data.md#44-consolidated-build-script) also go
  through a value channel and an `input: path`.

  **Worked example — samplesheet reads.** Raw reads are never handed to
  a process as path strings recovered from a node-local file. C1
  (`VALIDATE_SAMPLESHEET`) only checks the CSV and emits each `reads`
  entry verbatim; `main.nf` then resolves each entry on the launcher
  into a `Path` (`file(params.data_dir).resolve(rel)` or
  `file('az://…')`), checks every one exists in a single pass, and
  passes the `Path`s to `PARSE_SAMPLESHEET`. Nextflow stages only the
  listed files, and only into that process. `PARSE_SAMPLESHEET` must
  remain the **single consumer of raw reads**: everything downstream
  consumes its concatenated output, so no other process may take raw
  reads as input.
