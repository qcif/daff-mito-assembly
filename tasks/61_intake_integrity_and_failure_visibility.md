# Task 61 — Intake integrity check, and no sample disappears silently

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** open. The status-vocabulary decision in §2.2 is made
(`invalid_input`, 2026-10-14). Best done after task 60_ci_unit_tests.md,
so the new Python is gated by CI. **Task
70_pipefail_and_task_error_surfacing.md must land before, or with,
§2.3(a).** Without `pipefail`, a FILTLONG failure emits an empty FASTQ.
At the moment only the resulting NanoPlot crash turns that into
`error`, and once §2.3(a) makes NanoPlot optional the sample would
become a silent 0× coverage `fail`. §2.3(c)'s `onComplete` handler
should extend the one task 70 §2.2 may add, not duplicate it. Task
71_duplicate_read_id_dedup.md extends §2.1's intake check with
duplicate read-ID detection.

## Overview

On 2026-10-13 a client upload (Lomandra) exposed two linked failures.

**1. The input file was corrupt, and nothing checked it.** The
uploaded `barcode17_P2b_Walsh.fastq.gz` (5.4 GB) is valid for its
first 454,438 reads. After that the gzip stream breaks, and later
records have sequence and quality strings of different lengths.
Nothing checks the reads on the way in. `CHOPPER`, `FILTLONG`,
`RECRUIT` and `METAFLYE` all ran on whatever part of the file decoded.
Only `NANOPLOT_RAW` noticed, and it crashed.

**2. The crash made the sample vanish.** Since task 49, every
per-sample process runs with `errorStrategy = 'ignore'`
(`conf/base.config`). That is meant to stop one sample's failure
killing its siblings (CONSTITUTION principle 8). It works for the
siblings, but an ignored task emits nothing. `main.nf` assembles
COLLATE's input with a chain of `.join(..., by: 0)` calls, and a join
silently discards a sample that is missing from either side. The
crashed sample is joined against `NANOPLOT_RAW.out.reports` on both the
`ok` and the `fail` branch, so it never reached COLLATE.

C7 (`bin/run_report.py`) is designed for this case. It compares the
samplesheet with the `metadata.json` files it receives, and reports any
sample without one as `error`. But this was a one-sample run.
`COLLATE.out...collect()` on an empty channel emits nothing, so
`RUN_REPORT` never ran either. The run had no per-sample bundle and no
run report, and Nextflow said "completed" with exit code 0. That is the
worst kind of failure for a deployed prototype, because nobody would
notice. It breaks principle 7 (a degraded sample must never look like a
non-event) and principle 8 (a per-sample failure becomes a marker that
goes on to COLLATE).

This task:

- adds an integrity check at intake, so a corrupt or truncated file
  fails loudly at the first stage and costs no compute. The sample's
  report then says which file is bad and where.
- makes sure any sample that drops out for any reason is still reported:
  the run report is always produced, and the operator gets an end-of-run
  signal.
- stops a failure in QC-only stages (NanoPlot) from dropping the
  sample.

## 1. Constraints

- **Principle 7 / rule 18:** a corrupt input, a crashed stage and a
  sample that is genuinely negative must be distinguishable at every
  surface: `metadata.json`, `report.html`, `run-report.html` and
  `run_manifest.json`.
- **Principle 8:** the intake check's verdict is *data* (a status
  file, exit 0), not a Nextflow error. Siblings are unaffected.
- **Single consumer of raw reads** (spec §1a, worked example).
  `PARSE_SAMPLESHEET` is the only process that sees raw reads. The
  check runs there, on files already staged to that node. No new
  process may take raw reads.
- **Rule 14:** the check is `bin/` Python, stdlib only, with 100 %
  branch coverage. No new third-party dependency.
- **Rule 12:** `PARSE_SAMPLESHEET` and `VALIDATE_SAMPLESHEET` currently
  run in `python:3.12-slim` by tag, not digest (`conf/containers.config`,
  marked TODO). Moving `PARSE_SAMPLESHEET` to the digest-pinned
  `neoformit/daff-wf5-scripts` image is part of this task. Moving
  `VALIDATE_SAMPLESHEET` too is a one-line change, so do both.
- **Cost:** the check reads every byte of every input. Measure it
  (§3.3). It must cost no more than `NANOPLOT_RAW` already does on the
  same file.

## 2. Scope

### 2.1 Intake check (new `bin/intake_reads.py`, run in `PARSE_SAMPLESHEET`)

`PARSE_SAMPLESHEET` today is `cat input*/* > <id>.reads.fastq.gz`.
Replace it with one stdlib script that validates each input, writes the
combined output, and writes a status file.

Per input file, in samplesheet order:

1. **Compression.** Detect gzip by magic bytes, not by extension.
   Decompress the whole stream. A truncated or corrupt stream is an
   error, as is trailing garbage after the last complete member.
2. **Record structure.** Records are four lines. The header starts
   with `@` and the separator starts with `+`. `len(seq) == len(qual)`.
   The last record is complete. An empty file, or one with zero
   records, is an error.
3. **First error wins per file.** Stop reading a file at its first
   error. A 5 GB file should not be read to the end just to collect
   more errors. Still check every file, so a row with several bad
   files names all of them.

Record for each file: its name, whether it is compressed, the number of
records read, the number of bytes, and the error if any (its kind, the
1-based record number, and a short detail). Error kinds are
`gzip_stream`, `truncated_record`, `length_mismatch`, `bad_header`,
`bad_separator` and `empty`.

**Output**, only if every file passes:

- gzip inputs are byte-concatenated, as now. A multi-member gzip file
  is valid.
- **Plain `.fastq` inputs are compressed on the way through.** C1
  accepts `.fastq`/`.fq`, but `PARSE_SAMPLESHEET` names its output
  `.fastq.gz` regardless, and `CHOPPER` runs `zcat` on it. A plain
  input therefore probably fails at `CHOPPER` today. **Confirm this
  first** with a small plain-FASTQ stub, and record the result.

**Status file**, written whatever the outcome (`intake.json`):

```
{
  "$schema": "wf5/intake/v1",
  "status": "valid" | "invalid",
  "files": [
    {"name": "barcode17_P2b_Walsh.fastq.gz", "compressed": true,
     "records": 454438, "bytes": 5412345678,
     "error": {"kind": "gzip_stream", "record": 454439,
               "detail": "invalid compressed data after 1.98 GB"}}
  ],
  "total_records": 454438
}
```

The numbers above are illustrative. The script always exits 0 once it
has written `intake.json`. A non-zero exit means the script itself
broke (for example, a missing argument), which correctly stays a
Nextflow error.

**Performance.** Python's `gzip` decompresses in C. The cost is
line iteration. If §3.3 shows the stdlib reader is too slow, the
fallback is to read decompressed text from a `zcat` pipe, keeping the
record checks in Python. Do not reach for a third-party parser.

### 2.2 Routing an invalid sample, and what it is called

`PARSE_SAMPLESHEET` emits `intake.json` alongside the reads. The reads
output becomes optional. In `main.nf`, branch on `intake.json` the same
way `COVERAGE_GATE` output is branched: an explicit allowlist, with
anything unrecognised failing closed (spec §2.1.3, rule 18).

- `valid` → `NANOPLOT_RAW` and `CHOPPER`, as now.
- `invalid` → straight to `COLLATE` with a minimal bundle. No QC, no
  recruitment and no gate. COLLATE gains `intake.json` as an input on
  every branch, so the per-sample `metadata.json` always records
  intake. C6 classifies the sample from it before looking at
  anything else.

**Status vocabulary: a new value, `invalid_input`** (user decision,
2026-10-14). The operator action is "fix and re-upload the file". That
is neither `fail` (a coverage decision about a good sample) nor `error`
(spec §6a: "the pipeline breaking", which sends someone looking for a
bug). Reusing `error` was rejected because the operator could not tell
"your upload is broken" from "our pipeline is broken" without reading
the reason.

This touches `assets/sample_metadata.schema.json`,
`assets/run_manifest.schema.json`, C6 and C7 constants, the report
labels and colours, and spec §6a and §2.2 (C6/C7 rows). It is a
client-visible vocabulary change in `metadata.json`.

`sample_status_reason` names each bad file, the error kind and the
record number in plain words, for example "barcode17_P2b_Walsh.fastq.gz
is corrupt: the gzip stream breaks after read 454,438. Re-export or
re-upload the file." The per-sample `report.html` shows this in the
same style as the coverage-gate failure page.

### 2.3 No silent drop

Three changes. Each closes a different path to "a sample vanished and
the run exited 0".

**(a) QC reports do not gate COLLATE.** `NANOPLOT_RAW` and
`NANOPLOT_CLEAN` produce reports and nothing else. Join them into
`ch_ok_inputs` and `ch_failed_inputs` so that a missing report gives
an empty slot, not a lost sample. C6 already treats a missing NanoPlot
directory as absent (`parse_nanostats`). Make sure `metadata.json`
says the QC report is unavailable (null `read_qc` sections are not
enough on their own; add a short reason) and that `report.html` shows
it. Do not show an empty table.

> **Pitfall.** `join(..., remainder: true)` also emits items that
> appear only on the *right*-hand side. Every sample has a NanoPlot
> report, so on the `ok` branch the `fail`-branch samples' reports
> would leak through with a null left-hand side, and the other way
> round. Filter those out after the join, and cover it with the
> two-sample stub test in §3.2. The existing `BIN_TARGET.out.isoforms`
> remainder join avoids this only because isoforms exist only on the
> ok branch.

Audit every other input of the two joins. List each one in Outcomes
as *essential* (the bundle cannot be built without it) or *optional*.
Make optional inputs remainder joins only where C6 already copes with
their absence. Record the rest in tasks/todo.md. Don't widen this
task's scope to make every stage optional.

**(b) `RUN_REPORT` always runs.** Give the collected metadata channel
an empty-list fallback so `RUN_REPORT` runs with zero bundles. Every
samplesheet row then reconciles as `error` (or `invalid_input`). Check
that C7 and the `path` input both handle an empty list. Add a C7 unit
test for "zero bundles, two samplesheet rows".

**(c) The operator is told at the end of the run.**

- `main.nf` gets a `workflow.onComplete` handler. When
  `workflow.stats.ignoredCount > 0` it logs a prominent warning: how
  many tasks failed and were ignored, that the affected samples are
  listed in `run-report.html`, and the `nextflow log` command that
  lists failed tasks with their sample tags. It changes no exit codes.
- `deploy/azure/run-wf5-azure.sh`: after Nextflow exits 0, read the
  local `run_manifest.json` and print a per-status count. If any sample
  is `error` or `invalid_input`, print those sample ids and their
  reasons, and exit with a **distinct non-zero code** (for example 3,
  documented in the usage text). In `--integration` mode, run
  `assertions.sh` first so its output is not hidden.

C7's `error` reason stays generic. C7 cannot see which stage failed.
The `nextflow log` hint from (c) is the route to that.

## 3. Testing

### 3.1 Unit tests

New `scripts/tests/test_intake_reads.py`, at 100 % branch coverage.
Build small inputs in `tmp_path` with stdlib `gzip`:

- A valid gzip input, a valid plain input, and a valid multi-member
  gzip input all pass. Record counts are correct. The output decompresses
  to exactly the input records in input order. Gzip inputs are
  byte-identical in the output.
- A mixed row (one gzip and one plain) gives a valid gzip output.
- Each error kind is reported with the right record number: truncated
  gzip (cut the bytes mid-stream), corrupt gzip (flip a byte inside the
  compressed data), trailing garbage, a record with seq and qual of
  different lengths, a missing `@`, a missing `+`, a final record
  missing its quality line, an empty file, and a gzip file whose
  content is empty.
- With two inputs where the second is bad, the first is fully counted,
  the second has the error, the status is `invalid`, and no reads
  output is written.
- Gzip detection uses magic bytes. A gzip file named `.fastq` and a
  plain file named `.fastq.gz` are both handled correctly.
- Reading stops at the first error (assert the reader is not consumed
  past it, for example with a generator input that raises if pulled
  further).

C6 (`scripts/tests/test_collate.py`):

- `intake.json` invalid → `invalid_input`, a minimal bundle, and a
  reason naming the file and record.
- `intake.json` valid with NanoPlot directories missing → the normal
  status, and the QC section marked unavailable.
- Unknown intake status → fails closed, the same pattern as the
  existing unknown-gate-status test.

C7 (`scripts/tests/test_run_report.py`): zero bundles with two rows
gives two `error` samples and correct counts. `invalid_input` is
counted and labelled.

Report (`scripts/tests/test_report.py`): the intake-failure page and
the "QC unavailable" state render.

### 3.2 Stub run (CI)

The stub `VALIDATE_SAMPLESHEET` emits one hard-coded sample
(`STUB-01`). Add a second (`STUB-02`, a different `assembly_target`) so
the stub run covers two samples. Update the `tests.yml` output-layout
check to expect both.

Add a failure-injection stub test as a second step in the
`nextflow-test` job. It runs with an extra config file under
`tests/stub/` that makes `NANOPLOT_RAW` fail for `STUB-02` only.
Prefer a dynamic directive in that config (for example a
`beforeScript` closure on `meta.sample_id`) over a test hook in the
module. Fall back to a `task.ext` flag read by the stub block only if
the directive approach doesn't work, and say which was used. Assert:

- Both samples have `metadata.json` and `report.html`. `STUB-02`'s QC
  section is marked unavailable.
- `run_manifest.json` lists both samples.
- `nextflow` exits 0, and the log contains the §2.3(c) warning.

Add a second injection in which every sample fails at `NANOPLOT_RAW`
**and** at `CHOPPER`. Assert that `run_manifest.json` still exists and
lists every sample as `error`. This is the Lomandra shape.

### 3.3 Integration and real data

- **Integration fixtures:** the three fixtures are valid gzip files,
  so intake must report `valid` with the correct record count for
  each. Assert this in `tests/integration/assertions.sh` (record
  `intake.json` `total_records` per fixture in its
  `tests/integration/expected/*/` directory). Everything downstream
  stays unchanged. Run the full integration profile and confirm the
  existing assertions pass.
- **Corrupt-input fixture, generated, not uploaded:** have
  `tests/integration/fetch_fixtures.sh` make
  `INT-CORRUPT-01.fastq.gz` from the first N bytes of an existing
  fixture (a deterministic `head -c`), and add it to the integration
  samplesheet as a fourth row. It stops at intake, so it costs seconds.
  Assert `invalid_input`, that the reason names the file, and that
  siblings are unaffected.
- **Timing:** measure `PARSE_SAMPLESHEET` wall time against
  `NANOPLOT_RAW` on the largest fixture, and on a ≥5 GB file if one is
  available (ask the user for the Lomandra file's location). Record
  both in Outcomes.
- **Lomandra re-run** (if the user provides access): the corrupt file
  stops at intake with the reason above. Then on Azure (task 57's
  launcher) the launcher prints the per-status summary and exits with
  the distinct code.

## 4. Spec updates

- spec/01-pipeline-flow.md §1 diagram and spec/02-stages.md: stage 0
  now validates reads integrity. Add a C1b row (or extend C1) for
  `intake_reads.py` to §2.2, with its image and the status file.
- spec/02-stages.md §2.1.4: replace the stale "`errorStrategy 'ignore'`
  only on this process" wording with the task 49 default, and add the
  three §2.3 guarantees. Note that an ignored task's outputs are absent,
  so every join into COLLATE is either essential (sample becomes
  `error` through C7) or remainder (sample proceeds without it).
- spec/06a-reports.md: the status list (§6a run summary), the
  intake-failure page, and "QC unavailable".
- Add `invalid_input` to the status vocabulary
  wherever `error` is listed (spec §6a, §2.2 C6/C7 rows), and to both
  schemas.
- spec/04a-azure-blob-storage.md or the launcher docs (wherever
  `run-wf5-azure.sh` usage is described): the new exit code.

## 5. Out of scope

- Naming the failed stage in C7's `error` reason. That would need
  per-task failure markers or trace parsing.
- FASTQ content checks beyond structure (quality-encoding range, read
  length distribution). Duplicate read IDs are task 71.
- Making every downstream stage optional in the COLLATE join. §2.3(a)
  audits the joins and records the rest.
- Re-running the full Lomandra batch beyond the corrupt sample.

## 6. Checklist

- [ ] Plain `.fastq` + `CHOPPER` `zcat` behaviour confirmed and recorded
- [ ] §2.1 `intake_reads.py` + `PARSE_SAMPLESHEET` rewired; both
      samplesheet processes on the pinned scripts image
- [ ] §2.2 invalid branch → COLLATE; C6 and report handle it
- [ ] §2.3(a) NanoPlot remainder joins with right-only filter; join
      audit recorded
- [ ] §2.3(b) RUN_REPORT runs on zero bundles
- [ ] §2.3(c) `onComplete` warning; launcher summary and exit code
- [ ] §3.1 unit tests at 100 % branch; flake8 clean
- [ ] §3.2 two-sample stub plus both injections in CI
- [ ] §3.3 integration green; corrupt fixture asserted; timings recorded
- [ ] §4 spec updated
- [ ] Outcomes section written
