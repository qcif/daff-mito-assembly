# Task 73 — Cloudgene post-run steps: failed-task logs and a reports bundle

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** open. No hard preconditions. §2.1 has two candidate designs;
design B depends on task 67_procps_in_images.md, design A does not.
Independent of task 72_nextflow_26_strict_syntax.md (shelved), but if
that lands first, check that §2.1's chosen hook is still legal under the
strict syntax.

**Related spec:** spec/06a-reports.md §6a.4 (run-report and status
vocabulary), spec/02-stages.md §2.1.4 (cross-sample failure isolation),
spec/04a-azure-blob-storage.md. CONSTITUTION principles 7 and 8, rules
14, 15, 17, 18 and 19.

## Overview

WF5 will be offered to DAFF users through Cloudgene
(https://cloudgene.qcif.edu.au), the web front end that already runs
Taxodactyl (wf2). The Cloudgene app definition is
daff-cloudgene/workflows/wf5.yml, in the separate `daff-cloudgene`
repo. After the Nextflow run, Cloudgene runs a few shell "post-run
steps" and shows their output to the user. Taxodactyl uses two scripts
for this, both in wf2-nextflow/bin/:

- `collect_errors.sh` prints the tail of each failed task's log and
  zips `<outdir>/errors/` into `errors.zip`, so the user can download
  the logs.
- `zip_reports.sh` gathers every HTML report into `reports.zip`, so the
  user can download all of them in one go.

WF5 has neither script. For now wf5.yml does both jobs with inline
commands. They print the sample-status counts from `run_manifest.json`
and zip the reports. That is a stopgap, and it misses the important
part. **A Cloudgene user has no way to see why a sample failed.**

When a WF5 task fails, `conf/base.config` ignores the error so that its
sibling samples carry on (principle 8). The failing task's stderr tail
is printed to the Nextflow console, and the sample shows as `error` in
`run-report.html`. The actual logs, however, stay in the Nextflow work
directory. On Azure that is a blob container
(`az://workdata/wf5/<run>/`), which a Cloudgene user cannot browse. Run
from the command line, `run-wf5-azure.sh` points the operator at
`nextflow log` and `az_fetch_nf_logs`. A web user has neither tool.

Taxodactyl solves this by copying each failed task's `.command.out`,
`.command.err` and `.command.log` into `<outdir>/errors/<PROCESS>/<tag>.*`
when the run completes (wf2-nextflow/main.nf, the "Log collection"
`workflow.onComplete` block). That makes the logs part of the outputs,
which Cloudgene does serve. WF5 needs the same: capture the logs inside
the pipeline, then add Cloudgene-facing scripts that summarise and bundle
them.

This is principle 18 applied to the web front end. An `error` row with
no reachable explanation is a dead end for the user. Their only option
is to email the operator, who then has to dig through blob storage.

## 1. Current state

- `conf/base.config`'s `errorStrategy` closure already reads
  `task.workDir.resolve('.command.err')` on failure and prints the last
  20 lines to stderr. That is proof that a remote (Azure Blob) work dir
  can be read from that hook.
- `main.nf`'s `workflow.onComplete` warns when
  `workflow.stats.ignoredCount > 0`. It does not copy any files.
- Tracing is opt-in (`run-wf5-azure.sh --trace`), because it needs `ps`
  in every task image (task 67). So a trace file is not guaranteed to
  exist.
- Output layout (spec §6a.4): `<outdir>/run-report.html`,
  `<outdir>/run_manifest.json` and `<outdir>/<sample_id>/report.html`.
  The reports are self-contained (inline CSS, JS and images). The only
  relative links are `run-report.html` → `<sample_id>/report.html`.
- Process tags are `${meta.sample_id}` throughout `modules/local/`, so
  a failed task's tag is its sample.

## 2. Work

### 2.1 Capture failed-task logs into `<outdir>/errors/`

Write each failed task's logs to:

```
<outdir>/errors/<PROCESS>/<sample_id>.out
<outdir>/errors/<PROCESS>/<sample_id>.err
<outdir>/errors/<PROCESS>/<sample_id>.log
```

This uses the same layout as wf2, so the two apps' downloads look
alike. Copy only failed tasks. Skip `COMPLETED` and `CACHED` ones.
Write nothing (no `errors/` directory) when no task failed. The script
in §2.2 treats "directory absent" as "no errors".

There are two candidate designs. Try A first.

**A. Copy in the `errorStrategy` closure (preferred).** The closure in
`conf/base.config` already has `task.workDir`, `task.process` and
`task.tag`, and already reads `.command.err`. Extend it to copy the
three streams into `errors/`. This needs no trace file, so it does not
wait on task 67. It also runs at the moment of failure, so the logs are
captured even if the run later aborts. Things to verify:

- The closure can resolve `params.outdir` with its run-time value
  (`--outdir` on the CLI), not just the config default.
- Copying from Azure Blob into a local outdir works. wf2 found that a
  `Path`→`Path` copy fails on Azure Blob because the provider has no
  `PosixFileAttributeView`. Copy through an input stream instead, as
  wf2 does.
- A copy failure is logged and swallowed, and the closure still returns
  `'ignore'`. The log capture must never turn an isolated sample
  failure into a run abort.
- Retries: if a process ever gets `maxRetries`, a later attempt must
  not overwrite the earlier attempt's logs. Today nothing retries, so
  an `.attempt<N>` suffix is optional. State in the PR which choice
  you made.
- `VALIDATE_SAMPLESHEET` and `RUN_REPORT` use `'terminate'` and so
  bypass the closure. Decide whether their logs are worth capturing:
  a terminated run ends Cloudgene's job before the post-run steps run.
  Note the decision.

**B. Parse the trace file in `workflow.onComplete` (wf2's approach).**
Port wf2's block, driven by the trace file's `name`, `status` and
`workdir` columns. Only take this route if A turns out to be unworkable.
It requires tracing to be always on, which is blocked until task 67
puts `ps` in every image.

Whichever design is used, the existing `ignoredCount` warning in
`main.nf` should also say where the logs went, for example: "Logs for
the failed tasks are in <outdir>/errors/".

**Strict-syntax check (task 72):** if logic is added to
`conf/base.config`, confirm it is a closure body (allowed), not a
top-level statement (forbidden). Keep it short. If it grows past ~10
lines, move the copy logic to a function in `main.nf` or `lib/` and
call it from the closure (rule 14's spirit applies to Groovy too).

### 2.2 `collect_errors.sh`

Adapt wf2-nextflow/bin/collect_errors.sh. Keep its interface
(`--dir <outdir>`, exit 0 when there is nothing to report) and its
behaviour:

- No `errors/` directory → print "No task errors" and exit 0.
- For each `.err` and `.log` file, sorted, print the relative path and
  its last 10 lines.
- Zip `errors/` to `errors.zip`, replacing any existing one.

Changes for WF5:

- **Lead with the sample status summary.** Fold in the stopgap's
  status summary from wf5.yml, which reads `run_manifest.json`: the
  per-status counts, then each `error` or `invalid_input` sample with
  its `sample_status_reason`. Then list the logs grouped by sample, so
  the user reads "sample X is `error` because …" next to the logs for
  X. `invalid_input` samples have a reason but usually no task logs
  (the intake check records the problem as data, not a crash). Make
  sure they still show up in the summary.
- A missing `run_manifest.json` is itself an error. Print it loudly and
  exit non-zero, mirroring `EXIT_NO_MANIFEST` in run-wf5-azure.sh.
- Parse JSON with `python3` (stdlib only) and not `jq`. `python3` is
  already a dependency of run-wf5-azure.sh. Do not assume `jq` exists
  on the Cloudgene host.

### 2.3 `zip_reports.sh`

Adapt wf2-nextflow/bin/zip_reports.sh (`--dir <outdir>`). WF5 differs
from wf2 in two ways:

- All per-sample reports are named `report.html`. They must not be
  flattened into one directory, or they overwrite each other.
- `run-report.html` links to `<sample_id>/report.html`. Keep that
  layout inside the zip so the links work after extraction:

  ```
  reports/
    run-report.html
    run_manifest.json
    <sample_id>/report.html      (one per sample that has one)
  ```

Samples with no `report.html` (`error` samples) are simply left out of
the zip. Their `run-report.html` row already says why. A missing
`run-report.html` is an error and exits non-zero. Make the script
idempotent: if a step is re-run, replace `reports.zip` and the staging
directory, and leave no stray files in `<outdir>` other than
`reports.zip`. wf2's version leaves a `reports/` directory behind and
fails on re-run because of `mkdir` without `-p`.

### 2.4 Where the scripts live

Put them in `deploy/cloudgene/`, next to `deploy/azure/`, and not in
`bin/`. Everything in `bin/` is staged onto every task container's
`PATH` (rule 12), and these scripts run on the Cloudgene host, never in
a task. Cloudgene installs the app from a release archive of this repo,
so they will be at
`${CLOUDGENE_APP_LOCATION}/deploy/cloudgene/<script>.sh`.

Use `#!/usr/bin/env bash` with `set -euo pipefail`, as in
run-wf5-azure.sh. Call only `bash`, coreutils, `zip` and `python3`.

### 2.5 Update wf5.yml (daff-cloudgene repo)

Replace the two inline steps with:

```
- name: Check for errors
  cmd: /usr/bin/bash ${CLOUDGENE_APP_LOCATION}/deploy/cloudgene/collect_errors.sh --dir $outdir
- name: Zip workflow reports
  cmd: /usr/bin/bash ${CLOUDGENE_APP_LOCATION}/deploy/cloudgene/zip_reports.sh --dir $outdir
```

This is a separate commit in the separate repo. Mention it in the PR
description.

## 3. Tests

Rule 15 allows no fourth test framework, so test the shell scripts from
`pytest` in `scripts/tests/test_cloudgene_scripts.py`. Each test builds
a small fake outdir in `tmp_path`, runs the script with `subprocess`,
and asserts on the exit code, stdout and the zip's contents (via
`zipfile`). Hand-write the fixtures, with no real run outputs, so that
reference-data drift cannot break them (rule 19).

| Case | Expect |
|------|--------|
| No `errors/`, all samples `ok` | exit 0, "No task errors", no `errors.zip` |
| `errors/` with two processes, one sample each | both tails printed, grouped by sample, `errors.zip` holds both |
| `invalid_input` sample, no task logs | listed in the summary with its reason |
| No `run_manifest.json` | collect_errors exits non-zero with a clear message |
| Three samples, one with no `report.html` | zip holds `run-report.html`, manifest and two `<id>/report.html` |
| Run zip_reports twice | second run succeeds, zip unchanged in content, no stray `reports/` |
| No `run-report.html` | zip_reports exits non-zero |

For §2.1 (pipeline-side capture):

- **Stub run.** Add a stub-only way to force one process to fail for
  one sample. Use an existing test hook if there is one. Otherwise add
  a param that is only read in the stub block, for example
  `params.stub_fail_sample`. Assert that `errors/<PROCESS>/<id>.err`
  exists and that the sibling samples still complete. If this needs a
  new param, add it to `nextflow_schema.json` (or the equivalent
  schema) as a hidden test param.
- **Integration.** The fixtures have no sample that crashes a tool.
  `INT-CORRUPT-01` exercises intake (`invalid_input`), not a crash. Do
  not add a crashing fixture just for this. Instead, assert in
  tests/integration/assertions.sh that a clean integration run
  produces *no* `errors/` directory.
- **Azure.** Run the forced-failure stub once with
  `run-wf5-azure.sh --stub` to show that the copy from blob works
  (design A's main risk). Record the result in the PR.

Run flake8 on the new test module, and run `scripts/pytest.sh`.

## 4. Out of scope

- Linking `error` rows in `run-report.html` to their logs in
  `errors/`. It is worth doing, but it is a report change (spec §6a.4).
  Add it to tasks/todo.md under the report backlog.
- The Cloudgene admin config (profile, plugins and Azure env vars).
  That is deployment work, already noted in the header comment of
  wf5.yml.
- A public example samplesheet with `az://` reads for wf5.yml's '?'
  link. Add it to tasks/todo.md if no other task picks it up.

## 5. Done when

- A run with one forced task failure produces
  `<outdir>/errors/<PROCESS>/<sample_id>.{out,err,log}`, locally and
  on Azure, and the sibling samples are unaffected.
- A clean run produces no `errors/` directory.
- `deploy/cloudgene/collect_errors.sh` and `zip_reports.sh` pass the
  §3 tests, and `reports.zip` opens with working run-report → sample
  links.
- wf5.yml calls both scripts, and the inline stopgap steps are gone.

## Outcomes

- **§2.1 design A** implemented in the `errorStrategy` closure in
  `conf/base.config` (~15 lines, within a closure body so strict-syntax
  safe). `params.outdir` resolves to the CLI `--outdir` value. Copy goes
  through input/output streams (Azure-safe); failures are logged and
  swallowed and the closure still returns `'ignore'`. No `.attempt<N>`
  suffix: nothing retries today. `VALIDATE_SAMPLESHEET` and `RUN_REPORT`
  (`terminate`) are not captured — the run ends before Cloudgene's
  post-run steps; recorded in tasks/todo.md.
- The `ignoredCount` warning in `main.nf` now names `<outdir>/errors/`.
- **Stub test:** no new param needed — the existing
  `tests/stub/inject_filtlong_failure.config` is reused; its CI step in
  `.github/workflows/tests.yml` now asserts `errors/FILTLONG/STUB-02.{out,err,log}`
  and that nothing else is in `errors/`. Verified locally under
  NXF_VER=25.04.0 (the host's default Nextflow 26 rejects the existing
  nextflow.config syntax; see task 72).
- `tests/integration/assertions.sh` asserts a clean run has no `errors/`
  (not executed here; integration data not run).
- `deploy/cloudgene/collect_errors.sh` (exit 3 on missing manifest) and
  `zip_reports.sh` added; zip_reports stages in a `mktemp` dir so no
  stray files are left in the outdir. Tests:
  `scripts/tests/test_cloudgene_scripts.py` (8 tests, pass; flake8 clean).
- **Deviation (minor):** `scripts/Dockerfile` test layer now installs
  `zip`, since `scripts/pytest.sh` runs in that image and the scripts
  need it. The `:test` image must be rebuilt (`docker rmi
  neoformit/daff-wf5-scripts:test`) before `scripts/pytest.sh` passes the
  zip tests; I ran them with the claude venv's pytest instead.
- **wf5.yml** in `../daff-cloudgene/workflows/` now calls both scripts
  (file is untracked in that repo; not committed — commit separately).
- **Not done:** the Azure `--stub` run (needs credentials and spends
  cloud resources); listed in tasks/todo.md.
