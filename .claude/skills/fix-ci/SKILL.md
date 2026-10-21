---
name: fix-ci
description: Diagnose and fix failing GitHub Actions workflows for this repo
---

Diagnose and fix failing GitHub Actions runs. "$1" may name a workflow or
run ID; otherwise find out what is red.

## Facts (verified; re-check if something fails)

- Workflows: `.github/workflows/tests.yml` (lint, python-tests,
  nextflow-test), `integration.yml` (real-data, ~40 min, nightly + on
  push to pipeline paths), `images.yml` (Docker build).
- **The user's "integration tests failing" may really be `Tests`.** At
  the time of writing, `Integration tests` was green and `Tests` red.
  Always check which workflow is failing rather than trusting the label.
- CI pins Nextflow **25.04.0**. The local `nextflow` is newer (26.x) and
  fails to parse `main.nf` (`Unexpected input: '+'`). Always prefix local
  runs with `NXF_VER=25.04.0`.
- `nextflow-test` has `needs: [lint]`. When lint fails, the stub job is
  *skipped*, not passed, so it may hide a second failure.

## Procedure

1. **Triage.** `gh run list --limit 10`, then `gh run view <id>` to find
   the failing job and step. Note any skipped jobs.
2. **Get the error.** `gh run view <id> --log-failed`. This has returned
   empty output before. If so, reproduce locally: extract the step's
   inline Python heredoc from `tests.yml` into the scratchpad and run it
   with `/home/cameron/.local/envs/claude/bin/python`, e.g.

   ```
   sed -n "/<step name>/,/^          EOF/p" .github/workflows/tests.yml \
     | sed -n '/<<.EOF./,/^          EOF/p' | sed '1d;$d' \
     | sed 's/^          //' > $SCRATCH/check.py
   ```

3. **Fix** using the catalogue below.
4. **Verify locally** by running every step of the failing job, plus any
   jobs it gates (`needs:`). For the stub job:
   `NXF_VER=25.04.0 nextflow run . -profile stub -stub-run --outdir tests/output`
   then confirm `STUB-01/{metadata.json,report.html}`,
   `run_manifest.json`, `run-report.html` exist. Clean up afterwards with
   `.claude/scripts/clean_nextflow_run.sh`.
5. **Push and confirm.** Only commit/push when the user approves, then
   `gh run watch` the new run until every job is green (including ones
   that were previously skipped).
6. **Extend the catalogue** below when you fix a new class of failure.

## Failure catalogue

- **Bare-params file check** (`Nextflow lint`): a `${params.x}` appears in
  a `script:` block of `modules/local/*.nf`. If `x` is a scalar
  (number/string), add it to `SCALAR_PARAMS` in `tests.yml`. If it is a
  path, wrap it: `${file(params.x)}`. Last seen: `partial_recovery_min_fraction`.
- **Container-coverage check**: a new `modules/local/*.nf` needs a
  `withName` entry in `conf/containers.config`.
- **No-latest-tag check**: pin a version tag in `conf/containers.config`.
- **flake8**: run `/home/cameron/.local/envs/claude/bin/flake8 scripts/ bin/ --exclude=scripts/report/static/`.
- **Stub end-to-end**: every process `stub:` block must create all
  declared outputs; check the layout assertions in `tests.yml`.
- **SIGPIPE under `pipefail`** (`Integration tests`, fixtures or
  Assertions step; log shows `Broken pipe` / `sort: write error`,
  exit 1/2/141): a pipeline whose last stage is `head` exits early.
  Use `head -n N < <(cmd)` or do the selection in awk. Seen in
  `fetch_fixtures.sh` (`zcat | head`) and `assertions.sh`
  (`sort | head -1`). `gh run view --log-failed` was empty; use
  `gh api repos/<owner>/<repo>/actions/jobs/<job_id>/logs`.
- **Integration tests**: `gh run download <id> -n integration-nextflow-diagnostics`
  (and `-n integration-outputs`), read `.nextflow.log` and
  `work/**/.command.err`, and compare against
  `tests/integration/assertions.sh`. Full reruns are too heavy for a
  workstation; use the `remote-run` skill.
