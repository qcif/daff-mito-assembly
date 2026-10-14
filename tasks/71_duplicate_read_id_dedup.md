# Task 71 — Detect duplicate read IDs at intake and deduplicate before QC

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** open. **Preconditions:** task
completed/70_pipefail_and_task_error_surfacing.md, so that a duplicate this task
misses still fails loudly, and task
61_intake_integrity_and_failure_visibility.md, whose
`bin/intake_reads.py`, `intake.json` and `invalid_input` status this
task extends. Design decisions in §2 were made by the user on
2026-10-14.

**Related spec:** spec/01-pipeline-flow.md §1 and §1a (single consumer
of raw reads), spec/02-stages.md §2 stages 0–3 and the C1/C6/C7 rows,
spec/06a-reports.md, spec/00-overview.md §0 (bundle layout).
CONSTITUTION principles 7 and 8, rules 12, 14, 16 and 18.

## Overview

Some client uploads contain the same reads more than once. Round 2's
`barcode21_P2_A_Tinaroo.fastq.gz` is the whole file written twice:
74,006 records, 37,003 distinct reads, and the two halves of the
decompressed file are byte-identical. This most likely comes from
concatenating a directory of FASTQs into a file inside the same
directory, or from merging the same `fastq_pass` folder twice.

Duplicates do two kinds of harm:

- **Filtlong refuses the file** (`Error: duplicate read name`). After
  task 70 that fails the sample as `error`, which is loud and correct,
  but the client gets no result from a sample whose reads are
  perfectly usable.
- **Where a tool *doesn't* refuse, the result is quietly wrong.**
  Every read counts twice. The coverage gate's estimate (C2) doubles,
  so a sample at 6× looks like 12× and passes the 10× hard floor. The
  assembler sees each read twice, as if it were independent
  support.

This task finds duplicates at intake, where every record is already
read, and repairs the only safe kind. The reads are reduced to one
copy per ID before any QC filtering, and the sample then runs
normally. The repair is recorded in `metadata.json` and shown in both
reports, so an auditor can see the input was altered (rule 16).

A duplicate ID whose copies are *different* reads can't be repaired
safely, because there is no way to tell which copy is real. That is a
broken upload, and it is reported as `invalid_input` (task 61).

### Design decisions (user, 2026-10-14)

- **Detect at intake, not after a FILTLONG failure.** An earlier idea
  was to catch Filtlong's error, run a dedup step, then retry
  FILTLONG. That was rejected. It depends on matching a tool's error
  text, adds a retry branch to the DAG, and leaves RECRUIT and the
  coverage gate exposed to duplicates whenever a stage order changes.
  Detecting at intake is deterministic and protects every stage.
- **Deduplicate by read ID with `seqkit rmdup -n`**, not bbtools
  `dedupe.sh`. `dedupe.sh` deduplicates by *sequence*. By default it
  also removes reads contained within longer reads, which is real data
  loss at organelle depth. And it would not fix a same-ID,
  different-sequence collision, which would still crash Filtlong.

## 1. Constraints

- **Single consumer of raw reads** (spec §1a). Detection lives in
  `PARSE_SAMPLESHEET` (task 61's `bin/intake_reads.py`). The new
  `DEDUPLICATE` process consumes `PARSE_SAMPLESHEET`'s output, not
  raw reads.
- **Principle 8.** A duplicate verdict is data in `intake.json`, not a
  Nextflow error.
- **Rule 18.** Only an *identical* duplicate is repaired. Anything
  ambiguous fails closed as `invalid_input`. The repair is never
  silent.
- **Rules 12 and 14.** `DEDUPLICATE` runs in a pinned seqkit
  biocontainer. All decision logic stays in `bin/intake_reads.py`,
  stdlib only, at 100 % branch coverage.
- **Cost.** The ID set adds memory to intake. Measure the peak RSS on
  the largest available file (§3.3) and confirm it fits
  `PARSE_SAMPLESHEET`'s resource label.

## 2. Scope

### 2.1 Detection (`bin/intake_reads.py`, extends task 61 §2.1)

While validating records, track read IDs across **every input file of
the sample**, not per file. A sample uploaded as two copies of one
file is the same defect.

- **Read ID** = the header's first whitespace-delimited token, minus
  the `@`. Confirm this matches how Filtlong, seqtk and minimap2 name
  reads, and record the check in Outcomes. If any of them uses the
  whole header, a duplicate under that rule must be caught too.
- For each ID, keep a compact fingerprint of the record (for example a
  short digest of sequence + quality). Store the fingerprint, not the
  sequence, so memory scales with read count rather than bases.
- Classify each repeat of an ID as:
  - **identical**: same fingerprint as the first occurrence;
  - **conflicting**: same ID, different fingerprint.

`intake.json` (schema `wf5/intake/v1`, bump the version if task 61 has
already shipped it) gains a sample-level `duplicates` block, shown
here with illustrative values:

```
"duplicates": {
  "identical": 37003,
  "conflicting": 0,
  "distinct_ids": 37003,
  "examples": ["167b535b-f9c5-4e93-8911-c7a50d670ddb", ...]
}
```

Cap `examples` at a small fixed number (constant, e.g. 5).

**Status**, after task 61's structural checks pass:

- no duplicates → `valid`, as now;
- identical duplicates only → `valid`, with `dedup_required: true`;
- any conflicting duplicate → `invalid`. The sample becomes
  `invalid_input`, and the reason names the file(s), the conflicting
  count and one example ID, in plain words. For example: "barcode21….
  fastq.gz contains 12 read IDs that appear more than once with
  different sequences, so the file cannot be repaired automatically.
  Re-export the reads from the sequencing run."

A structural error still wins. Detection only reports duplicates for a
structurally valid sample.

### 2.2 `DEDUPLICATE` process (new, stage 0b)

- Input: `PARSE_SAMPLESHEET`'s combined reads, for samples with
  `dedup_required: true` only.
- Tool: `seqkit rmdup -n` (keep the first occurrence of each ID),
  output gzipped. Use a pinned seqkit biocontainer. Prefer it over the
  bespoke coverage-gate image (rule 12).
- Emits the deduplicated reads plus a small stats file (records in,
  records out, removed).
- **Consistency check.** `removed` must equal `intake.json`'s
  `identical` count. A mismatch means detection and repair disagree
  about what a duplicate is. That must not proceed: fail the process
  so the sample is `error` (it is a pipeline defect, not a bad
  upload). Ten lines of shell or fewer is fine. Anything more goes in
  `bin/` (rule 14).

### 2.3 Wiring (`main.nf`)

Branch on `intake.json` with an explicit allowlist, failing closed
(spec §2.1.3, rule 18), as task 61 §2.2 does:

- `valid`, no dedup → `NANOPLOT_RAW` and `CHOPPER` as now;
- `valid`, `dedup_required` → `DEDUPLICATE` → `CHOPPER`.
  `NANOPLOT_RAW` stays on the **submitted** reads, because it
  describes what the client sent, duplicates included. The report
  says so (§2.4).
- `invalid` → task 61's minimal-bundle path to COLLATE.

Mix the two `CHOPPER` inputs back into one channel. Nothing downstream
of `CHOPPER` changes. `DEDUPLICATE`'s stats file joins into COLLATE as
an *optional* input (task 61 §2.3(a) pattern, right-only filter
included), because most samples have none.

### 2.4 Recording and reporting (C6, C7, reports)

- **C6 / `metadata.json`:** an `intake.duplicates` section with the
  counts above and whether deduplication ran, plus `DEDUPLICATE`'s
  in/out counts. The seqkit version goes in the tool-versions block
  (rule 16). Update `assets/sample_metadata.schema.json`.
- **`sample_status` is unchanged** by deduplication. A repaired sample
  is judged on its deduplicated reads like any other. The repair is a
  provenance fact, not a degraded result.
- **Per-sample `report.html`:** a visible notice in the QC section,
  for example "37,003 duplicate reads (every read appeared twice) were
  removed before analysis. Raw QC below describes the file as
  submitted." Same notice style as the existing warnings. Never shown
  when no duplicates were found.
- **`run-report.html` / `run_manifest.json`:** a per-sample flag that
  input was deduplicated, so a batch-level reader sees it without
  opening each sample. Update `assets/run_manifest.schema.json`.

## 3. Testing

### 3.1 Unit tests (100 % branch coverage)

`scripts/tests/test_intake_reads.py`:

- No duplicates → no dedup flag, zero counts.
- Whole-file self-concatenation (the barcode21 shape) →
  `identical == distinct_ids`, `dedup_required`.
- A duplicate spread across two input files of one sample → detected.
- Same ID, different sequence → `conflicting`, status `invalid`.
- Same ID, same sequence, different quality string → `conflicting`
  (the records differ, and nothing should guess which is right).
- Headers that differ only after the first whitespace token → treated
  per the rule chosen in §2.1.
- A structural error in a file that also has duplicates → the
  structural error is reported, and duplicates are not.
- `examples` is capped.

`scripts/tests/test_collate.py`: a sample that was deduplicated
records the section and keeps its normal status; a sample that
wasn't omits it (or records zeros, whichever the schema says); a
conflicting-duplicate `invalid` intake gives `invalid_input` with
the reason text.

`scripts/tests/test_run_report.py` and `scripts/tests/test_report.py`:
the dedup flag and the notice render, and are absent when there are
no duplicates.

### 3.2 Stub run (CI)

Extend the stub so that one stub sample's `intake.json` sets
`dedup_required`. Assert that `DEDUPLICATE` ran for that sample only,
that its output reached `CHOPPER`, and that its `metadata.json`
carries the duplicates section. Reuse task 61's two-sample stub.

### 3.3 Integration and real data

- **Generated fixture, not uploaded** (constitution rule 19: no new
  test data to drift). In `tests/integration/fetch_fixtures.sh`, build
  `INT-DUPID-01.fastq.gz` deterministically from the first N records
  of an existing fixture, followed by the same N records again. Choose
  N so the sample stops at the coverage gate, which costs seconds, not
  an assembly. Add it as a samplesheet row. Assert in
  `tests/integration/assertions.sh` that:
  - intake reports `identical == N`;
  - `DEDUPLICATE` emitted N records;
  - `FILTLONG` succeeded;
  - the sample's status is the gate's verdict, not `error`;
  - siblings are unaffected.
- **Real data:** re-run the round 2 plant samples on `daff-admin`
  (`/remote-run`, `tests/manual/samples_round_2_plant.csv`, fresh
  `--outdir`). Expect `barcode21_P2_A_Tinaroo` to deduplicate
  74,006 → 37,003 records and to reach a real verdict. Record its
  status, gate coverage and loci in Outcomes. The other four samples
  should be unchanged (no dedup notice).
- **Cost:** record intake's peak RSS and wall time with duplicate
  tracking on, against task 61's recorded timings, on the largest file
  available.

## 4. Spec updates

- spec/01-pipeline-flow.md §1 diagram: stage 0b `DEDUPLICATE`
  (conditional) between `PARSE_SAMPLESHEET` and `CHOPPER`. Note that
  `NANOPLOT_RAW` reads the submitted reads.
- spec/02-stages.md §2 table: a `DEDUPLICATE` row. Extend the C1b (or
  C1) intake row with duplicate detection and the
  identical/conflicting rule.
- spec/06a-reports.md: the dedup notice (per-sample) and flag
  (run-level).
- spec/00-overview.md §0: no samplesheet change. If the bundle layout
  gains a `qc/dedup/` intermediate, list it.

## 5. Out of scope

- Sequence-level deduplication (identical reads under different IDs).
  That needs a different tool and has a real false-positive cost at
  organelle depth.
- Telling the client *how* their file came to be duplicated. The
  report states what was found and what was done.
- Any retry of `FILTLONG` after failure (rejected, see Design
  decisions).

## 6. Checklist

- [ ] §2.1 read-ID rule confirmed against Filtlong/seqtk/minimap2;
      detection in `intake_reads.py`
- [ ] §2.2 `DEDUPLICATE` with pinned seqkit biocontainer and
      consistency check
- [ ] §2.3 wiring: allowlist branch, mixed `CHOPPER` input, optional
      COLLATE join
- [ ] §2.4 C6/C7, both schemas, both reports
- [ ] §3.1 unit tests at 100 % branch; flake8 clean
- [ ] §3.2 stub covers the dedup path
- [ ] §3.3 `INT-DUPID-01` asserted; barcode21 re-run and timings
      recorded
- [ ] §4 spec updated
- [ ] Outcomes section written
