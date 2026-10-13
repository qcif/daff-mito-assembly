# Task 70 — A failed tool inside a pipe must fail its process, and say why

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** open. No preconditions. **Do this before task
61_intake_integrity_and_failure_visibility.md §2.3(a)** (see §1.3).
Task 71_duplicate_read_id_dedup.md builds on it.

**Related spec:** spec/02-stages.md §2 stages 2–3 and §2.1.4,
spec/06a-reports.md (status vocabulary), spec/05-test-data.md §5a/§5b.
CONSTITUTION principles 7 and 8, rules 15 and 18.

## Overview

On 2026-10-14 a remote run of the round 2 plant samples (`/remote-run`,
outputs in `tests/manual/output/samples/`) reported
`barcode21_P2_A_Tinaroo` as `error` with the reason "no metadata.json
was produced". That was the right status, but it was reached by luck,
and the real reason was never shown anywhere.

What actually happened:

1. The client file contains every read twice (§1.1). Filtlong
   refuses input with duplicate read names. It scored about 36,500
   reads, printed `Error: duplicate read name: 167b535b-…` and
   exited 1, having written nothing.
2. `FILTLONG`'s script is `filtlong … | gzip > out.fastq.gz`. Nextflow
   runs process scripts with `bash -ue`, **without `pipefail`**, so
   a pipe's exit status is the exit status of its *last* command.
   `gzip` compressed an empty stream perfectly well, so the task
   exited 0 and passed a 20-byte empty FASTQ downstream.
3. `RECRUIT` and `COVERAGE_GATE` ran on zero reads. The gate wrote
   `status: fail` at 0.0× coverage. That would have shipped as a
   coverage-gate failure: "your sample had too little organelle DNA".
   That is the wrong finding, and it sends the client to resequence a
   sample that is fine.
4. It didn't ship, only because `NANOPLOT_CLEAN` happened to crash on
   the empty file ("No reads found in input"). Its output was then
   missing, the `main.nf` join into COLLATE dropped the sample, and C7
   (`bin/run_report.py`) reconciled it as `error` because it had no
   `metadata.json`.
5. The console only said ``Process `NANOPLOT_CLEAN
   (barcode21_P2_A_Tinaroo)` terminated with an error exit status
   (1) -- Error is ignored``. That names the wrong stage, gives no
   message, and Filtlong's real error appears nowhere in Nextflow's
   output.

This task fixes the general defect: a tool failing inside a shell pipe
is invisible. When `FILTLONG` (or any tool) fails, its process must
fail with the tool's own error. The sample is then reported as `error`
because the failure was *caught*, not because a later stage tripped
over the leftovers. The tool's message must also reach the operator on
Nextflow's console output.

### 1.1 The barcode21 input

`barcode21_P2_A_Tinaroo.fastq.gz` holds 74,006 records, which are
37,003 distinct reads, each appearing exactly twice. The first and
second halves of the decompressed file are byte-identical (same md5),
so the file is a self-concatenation. The other four plant files in
the same delivery have no duplicate IDs. Handling such files is task
71's job. This task only makes sure the failure is loud and correct.

### 1.2 Why this is not only a FILTLONG problem

There is no `pipefail` anywhere in the pipeline (no `process.shell`
setting in `nextflow.config` or `conf/`). Every module that pipes into
`gzip` has the same blind spot, for example:

- `CHOPPER`: `zcat | chopper | gzip`. A corrupt input or a chopper
  crash produces an empty or truncated output with exit 0.
- `RECRUIT`: `seqtk subseq | gzip`.
- `FILTLONG`: as above.

### 1.3 Interaction with task 61

Task 61 §2.3(a) makes NanoPlot outputs optional in the COLLATE joins,
so that a NanoPlot crash no longer drops a sample. That is right, but
step 4 above is the only thing that turned this case into `error`.
If task 61 lands first, this exact input becomes a coverage-gate
`fail` at 0.0×, which is a silent mis-signal (principle 7, rule 18).
Land this task first, or in the same change.

## 2. Scope

### 2.1 `pipefail` for every process

Set the process shell to bash with `-e -u -o pipefail` globally (the
`process.shell` directive in `nextflow.config`), so that any failing
command in a pipe fails the task.

Then audit every `script:` block in `modules/local/*.nf` for pipes
where a non-zero exit from an earlier stage is **expected**. Those now
fail spuriously. The classic cases:

- `… | head` (or any reader that stops early). The writer gets
  SIGPIPE and exits 141.
- `grep` with no match (exit 1) feeding a later command.
- Tools that return non-zero for a legitimate empty result.

List every pipe in Outcomes as either *safe under pipefail* or
*changed* (with the change). Don't suppress an exit status blindly
(`|| true`). Where a non-zero exit is genuinely benign, check that
specific status and say why in a comment.

Also audit the `bin/*.py` scripts that build shell or `subprocess`
pipelines (for example the `seqtk sample | gzip` subsample branch in
`bin/coverage_gate.py`). Python ignores the global shell directive, so
each stage's return code must be checked explicitly. Where one isn't,
fix it and cover it with a unit test that fakes a non-zero upstream
exit (the existing tool-boundary mocking pattern, spec §5a).

Confirm that every container the pipeline uses has `/bin/bash`. The
current default already relies on it, but check that nothing was
falling back to `sh`.

### 2.2 The tool's error reaches the operator

**Required behaviour.** When a per-sample task fails and its error is
ignored (`errorStrategy = 'ignore'`, `conf/base.config`), the
Nextflow console output (stdout or stderr of the `nextflow run`
process, not only `.nextflow.log`) must show:

- the process name and the sample tag,
- the exit status,
- the tail of the task's stderr (enough lines to include the tool's
  own message, for example `Error: duplicate read name: …`),
- the task's work directory.

Show this once per failed task, ideally as it happens. A single
end-of-run summary listing every failed task is also acceptable.
Exit codes stay as they are: the run still completes, and sibling
samples are unaffected (principle 8).

**Investigate first and record the findings.** Find out what Nextflow
25.04.0 prints for an ignored task (expected: only the `NOTE: … Error
is ignored` line), and what is available to each hook. Candidate
mechanisms, in order of preference:

1. An `errorStrategy` closure that reads the failed task's
   `.command.err` tail and logs it before returning `'ignore'`. Check
   which `task` properties the closure can see in 25.04.0 (for
   example `task.workDir`, `task.exitStatus`, `task.tag`).
2. A `workflow.onComplete` summary built from whatever Nextflow
   exposes about failed tasks. Task 61 §2.3(c) adds an `onComplete`
   warning: share the handler rather than adding a second one.
3. A documented `nextflow log <run> -f name,tag,exit,stderr`
   command printed at the end of the run. This is a fallback only,
   because it does not meet the requirement on its own.

Whatever is chosen must work under every profile, `docker` and the
Azure Batch executor included (rule 17). Reading a work dir's
`.command.err` from the launcher may not work on Azure without
staging it. If it doesn't, say so in Outcomes and note how the Azure
path differs.

### 2.3 Sample status

No change to the status vocabulary. A sample whose task failed and was
ignored produces no `metadata.json`, so C7 reports it as `error`,
which is correct: the pipeline could not process it.

Do **not** add per-stage failure detail to C7's reason. Task 61 §5
keeps that out of scope, and it still is. §2.2 is the route to the
detail.

## 3. Testing

### 3.1 Unit tests

Only where §2.1's `bin/` audit changes Python. Each fixed pipeline
gets a test that fakes an upstream non-zero exit and asserts the
script fails, rather than writing an empty output. Keep the existing
100 % branch coverage on every touched C-component.

### 3.2 Stub run (CI)

Add a failure-injection stub step, using the same mechanism task 61
§3.2 specifies (an extra config under `tests/stub/` with a dynamic
directive keyed on `meta.sample_id`). Make `FILTLONG` fail for one
sample in a way that models the real defect: the failing command is
**not** the last one in a pipe, and it writes a known message to
stderr. Assert that:

- `FILTLONG` is reported as failed (not `NANOPLOT_CLEAN` or a later
  stage) and no downstream process ran for that sample;
- the injected stderr message appears in Nextflow's console output;
- `run_manifest.json` lists the sample as `error` and every sibling
  is unaffected;
- `nextflow` exits 0.

If task 61's two-sample stub (`STUB-01`, `STUB-02`) is not in place
yet, add the second stub sample here and let task 61 reuse it.

A positive control proves the injection exercises `pipefail`: with
`pipefail` temporarily removed, the same injection should *not* fail
the task. Run it once by hand and record the result. Don't keep it in
CI.

### 3.3 Real tools

- Run filtlong from the pinned `FILTLONG` container on a three-record
  FASTQ with one duplicated name, through the module's own pipe.
  Record the exit status before and after the fix.
- Run the full integration profile. All existing assertions must
  pass. A spurious `pipefail` failure here means §2.1's audit missed
  something.
- Re-run the round 2 plant samples on `daff-admin` (`/remote-run` with
  `tests/manual/samples_round_2_plant.csv`, a fresh `--outdir`).
  Expect `barcode21_P2_A_Tinaroo` to fail at `FILTLONG`, to be
  `error` in `run_manifest.json`, and Filtlong's message to appear in
  `run_<name>.log`. The other four samples' statuses should be
  unchanged (`ok` for barcode22, `fail` for 18, 19 and 24).

## 4. Spec updates

- spec/02-stages.md §2.1.4: the shell runs with `pipefail`, and a
  tool failing anywhere in a pipe fails its process. State how an
  ignored task's error is surfaced (§2.2). Correct the stale
  "`errorStrategy 'ignore'` **only** to guard" wording if task 61
  hasn't already.
- spec/01-pipeline-flow.md §1a engineering constraints: add the
  `pipefail` rule, plus the rule that `bin/` scripts check every
  subprocess stage's return code.

## 5. Out of scope

- Deduplicating or otherwise repairing the input: task 71.
- Making NanoPlot optional in the COLLATE joins: task 61 §2.3(a).
- Naming the failed stage in `metadata.json` or the run report.

## 6. Checklist

- [ ] §2.1 `pipefail` set globally; every module pipe audited and
      listed in Outcomes
- [ ] §2.1 `bin/` subprocess pipelines audited; fixes unit-tested
- [ ] §2.2 mechanism chosen after investigation; ignored-task error
      text reaches the console; Azure behaviour recorded
- [ ] §3.2 FILTLONG failure-injection stub step in CI; manual
      no-`pipefail` positive control recorded
- [ ] §3.3 real-tool check, integration green, `daff-admin` re-run of
      round 2 plant recorded
- [ ] §4 spec updated
- [ ] flake8 clean
- [ ] Outcomes section written
