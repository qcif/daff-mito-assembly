# Task 49: round_2 insect remote validation run

## Overview

Every automated test surface this pipeline currently has runs against
*subsampled* fixtures: `-stub-run` touches empty files, and
`-profile integration` runs three small, pre-recruited SRA-derived
samples chosen to fit a 60-minute CI budget (spec §5.1). Neither surface
tells us how the pipeline behaves on a real, un-subsampled biosecurity
submission — full read depth, real basecalling artefacts, genuinely
unknown taxa. `reference-material/samples/round_2/insect/` is exactly
that: 9 real ONT FASTQ files from an actual round of client insect
samples, all unambiguously `animal_mt` by their directory placement.

This task runs the full production pipeline (real reference bundle, real
tool containers — not stubs, not the integration profile's small
fixtures) over all 9 files as one batch. Because the data volume and
compute needs exceed what should run on a developer workstation
unattended, the run happens on a remote host reachable via `ssh
daff-admin`, with every artefact the run produces or consumes — synced
refs, synced sample data, Nextflow's `-work-dir`, and the output
directory — confined under `/mnt/data/wf5/` on that host. Nothing this
task does should touch daff-admin outside that path.

This mirrors the precedent already in the repo at `tests/manual/`
(a two-sample real-data run against production refs, referenced in
task 39 §7) but at full batch scale and on a host sized for it, rather
than a developer laptop.

## Background

- CONSTITUTION hard constraint 2: `assembly_target` is a per-row,
  pre-assembly gate declared at intake. Every file under
  `reference-material/samples/round_2/insect/` gets `animal_mt` — the
  directory placement is the classification here, no further lookup
  needed.
- spec/00-overview.md §0 defines the samplesheet schema
  (`sample_id`, `assembly_target`, `reads`, plus optional
  `sample_info` / `sample_type` / `sample_receipt_date` /
  `storage_location`) and the `--data-dir`-relative path resolution
  rules this task's samplesheet must follow.
- spec/05-test-data.md §5.1 is the reason this run is worth doing at
  all: fixture numbers are a floor, not an estimate of real
  performance, precisely because the fixtures are small, old, and
  unreplicated. This run is real submission data at real depth and is
  the first chance to see whether thresholds calibrated on the
  fixtures (`coverage_limits`, `recruit_thresholds`,
  `bin_target_thresholds` — all flagged provisional at their
  definitions) hold up, at least for the `animal_mt` arm.
- `tests/manual/samples.csv` and `tests/manual/local_resources.config`
  are the existing pattern for a real-data run outside the fixture
  harness: a samplesheet naming real client files, plus a small
  process-resource override config layered on top of whichever
  container profile is in use. This task follows the same shape,
  scaled to 9 samples and a remote host.

## Scope

1. Build a samplesheet covering all 9 files in
   `reference-material/samples/round_2/insect/`, each with
   `assembly_target = animal_mt`.
2. Provision `/mnt/data/wf5/` on the `daff-admin` host: sync the pipeline
   checkout, the `refs/v2026.09_1` reference bundle, and the insect
   sample data + generated samplesheet there.
3. Execute the pipeline on that host, in a container profile (Docker or
   Singularity, whichever `daff-admin` supports), with `-work-dir` and
   `--outdir` both under `/mnt/data/wf5/`, pointed at the synced
   production refs.
4. Retrieve the run-level report and per-sample bundles back into this
   repo for review.
5. Record what the run showed — clean passes, `low_coverage` calls,
   `fail` / `no_assembly` / `no_barcode` calls, anything unexpected —
   against the three-signal negative-clarity model (CONSTITUTION
   principle 7).

## Out of scope

- `reference-material/samples/round_2/plant/` and the 14 loose
  top-level round_2 files. The plant subdirectory needs a hard-
  constraint-3 decision (one row per organelle wanted) and the loose
  files need kingdom confirmed against submission records before a
  samplesheet row can be written for them — neither is resolved by
  this task. A follow-up task can cover them once that classification
  work is done.
- No pipeline code changes are anticipated by this task. If the run
  surfaces a genuine defect, stop, document the failure mode and the
  sample(s) that hit it, and open a new numbered task to fix it rather
  than patching mid-run.
- No changes to `tests/integration/*` — this is not a CI fixture and
  should not become one without the fixture-drift consideration in
  CONSTITUTION rule 19 and the `tasks/todo.md` fixture-provenance item
  being weighed first (as task 39 §7 already flagged for one round_2-
  adjacent sample).
- Retuning any threshold in `nextflow.config` off this run's numbers.
  Per spec §5.1, tuning against a single batch is exactly the trap that
  cost task 28 §10.2 — this run is an observation, not a calibration
  exercise.

## Procedure (pseudocode / illustrative only — no code to write)

1. **Build the samplesheet.** Per the spec §0 schema, e.g.:

   ```
   sample_id,assembly_target,reads,sample_info,sample_type,sample_receipt_date,storage_location
   R2-INSECT-barcode02,animal_mt,barcode02_P2_A_AR21a-1-1-1.fastq.gz,...
   R2-INSECT-barcode03,animal_mt,barcode03_P2_A_MG250214-4-1.fastq.gz,...
   ... (one row per file in round_2/insect/, 9 total)
   ```

   `reads` paths are relative to whatever `--data-dir` ends up being on
   the remote host (mirrors `reference-material/samples/round_2/insect/`
   once synced).

2. **Provision the remote host.**
   - `ssh daff-admin`, confirm `/mnt/data/wf5/` exists (or create it),
     and check available disk space and CPU/RAM (`df -h`, `nproc`,
     `free -h`) against the refs bundle + insect data size plus working
     space for Flye assemblies.
   - Sync the pipeline checkout (`main.nf`, `modules/`, `bin/`, `conf/`,
     `assets/`, `nextflow.config`) to `/mnt/data/wf5/pipeline/`.
   - Sync `refs/v2026.09_1/` to `/mnt/data/wf5/refs/v2026.09_1/`.
   - Sync `reference-material/samples/round_2/insect/` and the
     generated `samples.csv` to `/mnt/data/wf5/data/round_2_insect/`.
   - If `daff-admin`'s CPU/RAM ceiling differs from the
     `process_medium` / `process_high` tiers in `conf/base.config`, add
     a small override config on the pattern of
     `tests/manual/local_resources.config` (sized from the values
     actually observed on the host, not guessed).

3. **Run.** From `/mnt/data/wf5/pipeline/` on `daff-admin`:

   ```
   nextflow run main.nf \
     --samplesheet /mnt/data/wf5/data/round_2_insect/samples.csv \
     --data-dir    /mnt/data/wf5/data/round_2_insect \
     --organelle_refs /mnt/data/wf5/refs/v2026.09_1/recruit \
     --protein_panel  /mnt/data/wf5/refs/v2026.09_1/proteins \
     --blast_db       /mnt/data/wf5/refs/v2026.09_1/validate \
     --annotate_refs  /mnt/data/wf5/refs/v2026.09_1/annotate \
     --refs_manifest  /mnt/data/wf5/refs/v2026.09_1/manifest.json \
     --outdir /mnt/data/wf5/output \
     -work-dir /mnt/data/wf5/work \
     -profile docker[,<remote resource override config if added>]
   ```

   Everything Nextflow reads or writes for this invocation must resolve
   under `/mnt/data/wf5/` — no default `-work-dir`, no default
   `--outdir` outside that path.

4. **Retrieve.** Copy `/mnt/data/wf5/output/` (per-sample bundles +
   `run-report.html`) back to `tests/manual/output/round_2_insect/` in
   this repo for review, following the existing `tests/manual/output/`
   convention.

5. **Review.** For each of the 9 samples, note its terminal status
   (`ok`, `low_coverage`, `fail`, `no_assembly`, `no_barcode`) and
   whether that status looks right given what's known about the
   sample. Cross-check the run-report.html summary against the
   per-sample reports. Flag anything surprising as a candidate
   follow-up task rather than resolving it here.

6. **Clean up.** Once outputs are retrieved and reviewed, clear
   `/mnt/data/wf5/work/` (per AGENTS.md's guidance on not letting stale
   `work/` data accumulate). Leave `/mnt/data/wf5/refs/` and
   `/mnt/data/wf5/data/` in place unless told to reclaim the space, and
   do not touch anything on `daff-admin` outside `/mnt/data/wf5/`.

## Acceptance criteria

- [x] All 9 files in `reference-material/samples/round_2/insect/` have
      a row in the samplesheet, `assembly_target = animal_mt`.
- [x] `/mnt/data/wf5/` on `daff-admin` holds the synced pipeline, refs,
      and sample data; no run artefact lands outside it.
- [x] The pipeline completes (or fails informatively) for every row in
      the samplesheet, using the real `v2026.09_1` reference bundle and
      real tool containers.
- [x] `run-report.html` plus all per-sample bundles are retrieved into
      `tests/manual/output/round_2_insect/` in this repo.
- [x] A short written summary exists (in the PR/commit description or a
      note alongside `tests/manual/output/round_2_insect/`) of
      per-sample terminal status and anything that warrants a follow-up
      task.
- [x] `/mnt/data/wf5/work/` is cleared after retrieval; no other change
      is made to `daff-admin` outside `/mnt/data/wf5/`.

## Outcomes

### Deviations from the brief

- **`--data-dir` in the brief's example invocation is wrong; the real
  flag is `--data_dir`** (`nextflow.config` declares `params.data_dir`,
  not a kebab-case alias). Used `--data_dir` throughout; the brief's
  §3 code block should be read with that correction.
- **Nextflow version.** The pre-installed `26.04.6` on `daff-admin`
  fails to parse `main.nf` at all (`Error main.nf:58:17: Unexpected
  input: '+'`, a stricter-parser rejection of an existing multi-line
  string concatenation) — reproduces identically on the local dev
  machine's `26.04.6`, so this is unrelated to the remote host or this
  task's changes. Installed Nextflow `25.04.0` instead (the version
  `.github/workflows/*.yml` already pins via `nf-core/setup-nextflow`),
  confined to `/mnt/data/wf5/nextflow`. Worth a follow-up note if
  anyone hits this again: the repo has no documented minimum beyond
  `nextflowVersion = '>=25.04.0'` in `nextflow.config`, which invites
  picking up whatever `get.nextflow.io` currently resolves to.
- **Docker data-root moved to `/mnt/data/docker` (host-wide change,
  outside `/mnt/data/wf5/`, done with explicit user approval).**
  `daff-admin`'s root disk (`/dev/vda2`, 30 GB) was already 93% full
  from unrelated prior work on this shared host (`neoformit/taxodactyl`
  images alone ~12 GB across tags) before this task touched anything,
  leaving no room to pull the real tool containers. Tried
  `-profile singularity` with its cache redirected under
  `/mnt/data/wf5/singularity_cache` first, to avoid any host-wide
  change — but the host's Singularity (`3.7.0`, ~2020) fails to
  convert at least one image's OCI manifest (`quay.io/biocontainers/
  nanoplot`: `FATAL: ... no descriptor found for reference ...`), a
  known old-Singularity/new-OCI incompatibility. With the user's
  explicit go-ahead, added `/etc/docker/daemon.json` (`{"data-root":
  "/mnt/data/docker"}`) and restarted the `docker` service; nothing
  under the previous root (`/var/lib/docker`) was deleted or migrated.
  Confirmed via `docker info` before proceeding.
- **`errorStrategy` gap fixed in the pipeline itself, not worked around
  per-run (explicit user direction, overriding the brief's default "no
  pipeline code changes" / "open a new task instead" guidance).** The
  first full-batch attempt aborted entirely after `R2-INSECT-barcode14`
  reached `METAFLYE`: Flye genuinely assembled 0 contigs and exited
  non-zero, and only `COVERAGE_GATE` had `errorStrategy 'ignore'` set —
  every other process defaulted to Nextflow's `terminate`, so one
  sample's real-tool failure killed the other 5 not-yet-started
  samples. Confirmed via `deploy/azure/run-wf5.sh` that production
  batches samples the same way (one samplesheet, one `nextflow run`),
  so this was a real production-relevant gap, not an artefact of this
  task's setup — and spec's `error` `sample_status` (`§6a.4`/`§2`,
  CONSTITUTION principle 7) already documents "a sample that crashed"
  as an expected, isolated outcome. Asked the user how to proceed;
  told to fix it as a real pipeline feature, generalised to every
  process rather than scoped to `METAFLYE`. Changed
  `conf/base.config`'s top-level `process {}` block to default
  `errorStrategy = 'ignore'`, with `VALIDATE_SAMPLESHEET` and
  `RUN_REPORT` (the two whole-run, non-per-sample gates) overridden
  back to `'terminate'`; removed the now-redundant per-process override
  from `modules/local/coverage_gate.nf`. Verified with a local
  `-profile stub -stub-run` before syncing to `daff-admin` and
  resuming. This is a real, intentional pipeline change from this
  task, included in the same commit.
- **`scripts/report/templates` / `scripts/report/static` were missed
  in the initial remote sync** (only `main.nf`, `modules/`, `bin/`,
  `conf/`, `assets/`, `subworkflows/` were synced per the brief's §2
  list, which itself omits `scripts/` even though `main.nf` sources
  `ch_report_templates`/`ch_report_static` from
  `${projectDir}/scripts/report/`). Every per-sample `report.html`
  rendered as a Jinja `TemplateNotFound` traceback page until this was
  caught and `scripts/` synced, then the batch re-resumed (cheap: only
  `COLLATE`+`RUN_REPORT` re-ran, everything else stayed cached). The
  brief's §2 sync list should be corrected to include `scripts/` for
  any future remote run.

### Per-sample results (9 files, all `animal_mt`)

| Sample | `sample_status` | Recruited coverage | Assembly | Barcodes recovered |
|---|---|---|---|---|
| barcode02 | `fail` | 0.13× | — (soft-failed pre-assembly) | — |
| barcode03 | `fail` | 5.05× | — | — |
| barcode04 | `fail` | 6.98× | — | — |
| barcode05 | `ok` | 143.35× | 1 contig | 6/6 loci |
| barcode06 | `low_coverage` | 26.18× | 1 contig | 4/6 loci |
| barcode08 | `no_assembly` | 35.57× (cleared the gate) | Flye produced 1 contig (3568 bp) but `BIN_TARGET` rejected it as off-target — a real, clean "recruited and assembled something, but not the target organelle" outcome, not a crash | — |
| barcode13 | `ok` | 1340.14× (subsampled to 297.84×) | 1 contig | 6/6 loci |
| barcode14 | *(crashed — no `metadata.json`)* | n/a | Flye genuinely produced 0 contigs (only 2.4% of reads aligned to each other) and exited non-zero | — |
| barcode18 | `fail` | 6.44× | — | — |

All three `fail` calls sit clearly under the 10× hard floor (0.13–6.98×)
— unambiguous, correctly-behaving coverage-gate soft-fails, not
borderline cases. `barcode06`'s `low_coverage` (26.18×, between the 10×
hard floor and 30× warn floor) is exactly the case that status exists
to describe, and still recovered 4/6 barcode loci. `barcode05` and
`barcode13` are clean high-confidence passes (143×, up to 1340× before
subsampling) with full 6/6 barcode recovery. `barcode08` and
`barcode14` are the two negative-clarity cases, and — per CONSTITUTION
principle 7 — they read as genuinely different things: `barcode08`
recruited and assembled real sequence that binning correctly
identified as not the declared target; `barcode14` never produced
assembleable sequence at all. Once `RUN_REPORT` (task 45) is real,
`barcode14` should surface as `sample_status: error` at the run level,
since it never produced a `metadata.json` — worth a spot-check when
that task lands.

### Known gap, not caused by this run

`run-report.html` and `run_manifest.json` are empty (0 bytes) in the
retrieved output — `RUN_REPORT`'s `script:` block is still a literal
`touch run_manifest.json run-report.html` stub regardless of profile
(`bin/run_report.py` doesn't exist yet). This is pre-existing and
already tracked under "Run-level provenance (task 45_run_report.md)"
in `tasks/todo.md`; not a new finding from this task, and not touched
here. Per-sample `report.html` (rendered by `COLLATE`, a separate
code path) is real and worked correctly once `scripts/` was synced.

### Follow-up

No new follow-up task filed. The one genuine pipeline defect this run
surfaced (`errorStrategy` gap) was fixed inline per the user's explicit
direction, and is included in this task's changes. The `RUN_REPORT`
stub gap was already tracked (task 45) before this run and needs no
new entry.
