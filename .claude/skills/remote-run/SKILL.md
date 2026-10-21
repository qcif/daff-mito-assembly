---
name: remote-run
description: Run the workflow on the daff-admin remote VM against real sample data
---

Run the pipeline on the `daff-admin` VM for "$1" — a real-data run at a
scale that shouldn't tie up a developer workstation.

**Everything this skill does stays under `/mnt/data/wf5/` on the remote.**
Never write outside it. The one standing exception is already in place
and needs no repeating (see Host facts).

## Host facts (verified; re-check rather than assume if something fails)

- Reachable as `ssh daff-admin` (4 vCPU, 7.8 GB RAM, `/mnt/data` ~980 GB).
- Nextflow lives at `/mnt/data/wf5/nextflow`, **pinned to 25.04.0** — not
  on `$PATH`, so always call it by that absolute path. Do not "upgrade"
  it: 26.x fails to parse `main.nf` outright (`Unexpected input: '+'` on
  an existing multi-line string concat). 25.04.0 is what
  `.github/workflows/*.yml` pins via `nf-core/setup-nextflow`.
- Docker's `data-root` is `/mnt/data/docker` (`/etc/docker/daemon.json`).
  The 30 GB root disk is shared with unrelated projects and stays near
  full — if you ever see `no space left on device` from a container
  pull, check `docker info | grep 'Docker Root Dir'` still points at
  `/mnt/data/docker` before doing anything else. Don't prune images:
  they belong to other projects on this host.
- `-profile singularity` does **not** work here (Singularity 3.7.0 can't
  convert current biocontainer OCI manifests). Use `-profile docker`.
- Reference bundles already staged under `/mnt/data/wf5/refs/`.
  Check `ls /mnt/data/wf5/refs/` for what's present before syncing 2 GB
  again.

## Procedure

1. **Samplesheet.** Build one per spec/00-overview.md §0 — required
   columns `sample_id`, `assembly_target`, `reads`; optional
   `sample_info`, `sample_type`, `sample_receipt_date`,
   `storage_location`. Unknown columns fail parse. `reads` paths are
   relative to `--data_dir`. Keep it in `tests/manual/` alongside the
   existing `samples*.csv`, and confirm `assembly_target` per row
   against the CONSTITUTION's hard constraint 2 (it's an intake
   declaration, never inferred later).

2. **Sync.** Run `bash .claude/scripts/sync_remote_pipeline.sh`. It
   mirrors every `${projectDir}`-resolved path and excludes
   `__pycache__`. Don't hand-roll the rsync — omitting `scripts/`
   (report templates) or `subworkflows/` is the classic failure, and a
   trailing slash on a directory arg flattens the layout. Add
   `--with-stub-fixtures` if you intend a remote stub run.

   It also copies `remote_resources.config` to `/mnt/data/wf5/` (the
   path step 4 passes to `-c`).

   Sample data and refs are separate (they're large and change rarely):
   ```
   rsync -a --progress <sample-dir>/ daff-admin:/mnt/data/wf5/data/<run-name>/
   rsync -a tests/manual/<sheet>.csv \
       daff-admin:/mnt/data/wf5/data/<run-name>/samples.csv
   rsync -az refs/<version>/ daff-admin:/mnt/data/wf5/refs/<version>/
   ```
   - **No `-z` for `*.fastq.gz`.** Re-compressing already-gzipped reads
     is pure CPU overhead — a 1.5 GB round took >10 min with `-z`.
     Run the transfer with `run_in_background` and wait on it rather than
     blocking on the 10-minute Bash timeout.
   - **Sanitise filenames on the way over.** Client files arrive with
     browser-download suffixes like `barcodeNN_X.fastq (1).gz`. Rsync
     those individually to a clean name
     (`... "<odd name>" daff-admin:.../barcodeNN_X.fastq.gz`) and use
     the clean name in the samplesheet; spot-check with `zcat | head -2`.
   - **Sample ID = the reads filename**, minus its FASTQ extension
     (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`): `barcode21_SAC1451-2-1.fastq.gz`
     → `barcode21_SAC1451-2-1`. Never shorten it to `barcodeNN`. Barcodes
     are reused across runs, so `barcode21` alone is ambiguous. Sanitise
     only as far as `sample_id` requires (`[A-Za-z0-9_.-]+`, spec/00 §0):
     drop browser suffixes like ` (1)` and replace any other disallowed
     character with `_`. For a multi-file sample, use the name the files
     share, and ask if they share none.

3. **Sanity-check the wiring first** when anything structural changed —
   cheap, and catches a bad sync before an hours-long run:
   ```
   ssh daff-admin "cd /mnt/data/wf5/pipeline && \
     /mnt/data/wf5/nextflow run main.nf -profile stub,docker -stub-run \
       --outdir /mnt/data/wf5/output_stub_check \
       -work-dir /mnt/data/wf5/work_stub_check"
   ```
   Then remove both `*_stub_check` dirs.

4. **Run detached.** The run outlives the ssh session, so `nohup` it and
   redirect to a log under `/mnt/data/wf5/` — don't hold it open in the
   foreground.

   **The flag is `--data_dir`, not `--data-dir`** (`nextflow.config`
   declares `params.data_dir`; the kebab form silently falls through to
   the profile default and fails on a path you never named).

   ```
   ssh daff-admin "cd /mnt/data/wf5/pipeline
   nohup /mnt/data/wf5/nextflow run main.nf \
     --samplesheet /mnt/data/wf5/data/<run-name>/samples.csv \
     --data_dir    /mnt/data/wf5/data/<run-name> \
     --organelle_refs /mnt/data/wf5/refs/<ver>/recruit \
     --protein_panel  /mnt/data/wf5/refs/<ver>/proteins \
     --blast_db       /mnt/data/wf5/refs/<ver>/validate \
     --annotate_refs  /mnt/data/wf5/refs/<ver>/annotate \
     --refs_manifest  /mnt/data/wf5/refs/<ver>/manifest.json \
     --outdir /mnt/data/wf5/output_<run-name> \
     -work-dir /mnt/data/wf5/work \
     -profile docker \
     -c /mnt/data/wf5/remote_resources.config \
     -with-report   /mnt/data/wf5/output_<run-name>/execution_report.html \
     -with-timeline /mnt/data/wf5/output_<run-name>/execution_timeline.html \
     -with-trace    /mnt/data/wf5/output_<run-name>/execution_trace.txt \
     > /mnt/data/wf5/run_<run-name>.log 2>&1 &
   disown"
   ```

   **Per-run `--outdir` and log.** Earlier rounds' results are left on
   the remote, so a shared `/mnt/data/wf5/output` would mix runs
   together. Check `ls /mnt/data/wf5/` first and never reuse another
   run's outdir. `\$` any shell variables you set inside the
   double-quoted ssh string.

   `remote_resources.config` caps `resourceLimits` to the host's real
   4 CPU / 6 GB — `conf/base.config`'s `process_medium`/`process_high`
   tiers (16/32 GB) exceed what this VM has.

5. **Wait on it properly.** Use a Monitor with an until-loop on the log
   rather than polling or sleeping:
   ```
   ssh daff-admin "until grep -qE 'Completed at|Execution cancelled|ERROR ~' \
     /mnt/data/wf5/run_<run-name>.log 2>/dev/null \
     || ! pgrep -f 'nextflow.*main.nf' >/dev/null; do sleep 30; done; \
     tail -c 3000 /mnt/data/wf5/run_<run-name>.log | tr '\r' '\n' | tail -40"
   ```
   Match failure signatures too, not just the success marker — silence
   must never be mistaken for progress. The `pgrep` arm catches a
   Nextflow that died without writing a marker (OOM-kill, host reboot).
   The log is ANSI progress output with `\r` redraws, hence the
   `tr`.

   A Monitor expires after 30 min max; a 4-sample round with one
   ~1 GB fastq takes ~40 min, so expect to re-arm at least once. On
   expiry, tail the log to confirm progress before re-arming.

6. **Retrieve and review.** Pull `/mnt/data/wf5/output_<run-name>/` into
   `tests/manual/output/<run-name>/` (gitignored, per existing
   convention). Report each sample's terminal status — `ok`,
   `low_coverage`, `no_assembly`, `no_barcode`, `fail`, `error` — and
   keep them distinct per CONSTITUTION principle 7. In particular
   `no_assembly` (assembled nothing, or binned off-target) and
   `no_barcode` (assembled fine, no locus extractable) are different
   findings, and `fail` (a deliberate coverage-gate decision) is not
   `error` (the pipeline broke). A sample that crashed outright leaves
   no `metadata.json` at all.

   Status lives in each `<sample>/metadata.json` as `sample_status` +
   `sample_status_reason`; also report `recovery.completeness` and its
   `reasons` alongside it (an `ok` sample can be a `partial` recovery);
   gate numbers are under `coverage.gate`,
   assembly summary under `assembly` (skip its `contigs` /
   `bin_metadata` — they're huge), loci under `barcodes.loci`. Print
   those fields selectively; dumping whole sections floods the output.

   `run-report.html` (batch summary, from `RUN_REPORT` / task 45) and
   `run_manifest.json` are real and are the first things to point the
   user at; per-sample `report.html` comes from `COLLATE`. A 0-byte run
   report now means something broke. The pull also brings back
   `execution_report.html` and `execution_timeline.html` (wall time,
   task table, resource use) and `execution_trace.txt`. Point the user at
   those too.

7. **Clean up.** `ssh daff-admin "rm -rf /mnt/data/wf5/work"` once
   outputs are retrieved, and trim stale `.nextflow.log.N` (keep
   `.nextflow.log` and `.1`/`.2`) in `/mnt/data/wf5/pipeline/`. Leave
   `refs/`, `data/` and `output_<run-name>/` in place unless told to
   reclaim the space.

## Resuming

`-resume` is cheap and correct for re-running after a report-side or
config fix — assembly stays cached. But per tasks/todo.md it may not
re-run a process after an in-place edit to `bin/*.py` (unconfirmed
whether script contents hash into the task hash), so don't trust it for
Python changes; delete the relevant work dir or run clean.

## If a sample fails

Per-sample crash isolation is on by default (`errorStrategy = 'ignore'`
in `conf/base.config`), so one sample's tool failure no longer aborts
the batch — `VALIDATE_SAMPLESHEET` and `RUN_REPORT` deliberately stay
fail-loud. A single ignored process is usually a real biological
outcome, not a bug: check whether the tool genuinely produced nothing
(e.g. Flye assembling 0 contigs) before treating it as a defect. Don't
patch pipeline code mid-run — document the failure and raise a task.

A task failing with "Command 'ps' required by nextflow to collect task
metrics cannot be found" means an image has lost `procps`. That is a
regression (task 67), not a sample failure.
