# Task 35 — `COVERAGE_GATE`: two floors, new status vocabulary, consumer migration

**Phase:** P4a (a revision to the P1 stage 6 gate, landing mid-P4a).
**Goal:** Replace `COVERAGE_GATE`'s single coverage floor with a **warn
floor** and a **hard floor**, so that under-covered samples are assembled
and flagged rather than discarded; rename the sample status vocabulary
accordingly; and migrate every consumer of the old vocabulary.

**Status:** implementation. The design decision is settled and already
recorded in
[spec §2.1](../spec/02-stages.md#21-coverage-gate-2-stage-6) — this task
implements the spec as written, it does not re-open it.

**Related:** revises the gate built in task 15 and corrected in task 25.
Triggers task 36 (`plant_mt` integration coverage), which depends on this
landing first. Partially supersedes the `INT-PLANT-01-mt` entry under
"Integration fixtures" in `tasks/todo.md`.

---

## 1. Overview

The coverage gate sits between read recruitment and assembly. It estimates
how deeply the target organelle was sequenced, and decides whether assembly
is worth attempting. Until now that decision was binary against a single
floor of 30×: above it the sample assembled, below it the sample was
"soft-failed" — marked as low coverage and skipped entirely, never reaching
the assembler.

That binary is too blunt for biosecurity work. A sample at 25× is not
hopeless; it is *degraded*. Running one such sample through the pipeline
with the floor temporarily lowered (a client `animal_mt` intercept,
`barcode06`, estimated 26×) produced a fragmented, non-circular assembly
covering about a third of the mitogenome — but it still yielded four
usable barcode loci including an intact `COX1`. Four loci, or even one, are
enough for Taxodactyl to make an approximate taxonomic assignment. Under
the single-floor gate that sample produced nothing at all, and the reads
were thrown away.

So the gate gains a second floor and the low end splits in three:

- Below the **hard floor** (10×) there is genuinely not enough read overlap
  for the assembler to build a graph. Skip assembly — the compute is
  wasted. This is the only terminal outcome.
- Between the hard floor and the **warn floor** (30×) assemble anyway, but
  mark the sample so nobody mistakes a partial recovery for a complete one.
- At or above the warn floor, business as usual.

Alongside this the three sample statuses are renamed to say plainly what
they mean: `fail` (gate rejected it), `low_coverage` (assembled, but under
the warn floor), `ok` (assembled normally).

**The trap in this task is the rename, not the logic.** `low_coverage`
already exists as a status and currently means *"rejected, no assembly
attempted"*. After this change the same string means *"assembled, passed
with a warning"* — very nearly the opposite. Nothing crashes when a stale
reader encounters it; it just quietly reports a successful assembly as a
failure, or vice versa. Section 6 exists because of this, and it is the
part of the task most likely to be under-done.

**Why the numbers are what they are.** 30× has one supporting observation
(the `barcode05`/`barcode06` pair in
[spec §2.1.6](../spec/02-stages.md#216-provenance-of-the-two-floors)); 10×
has none and is set from first principles. Both are explicitly provisional
and are swept properly in
[spec §9 item 2](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking).
This task implements the mechanism and the placeholder values. **It does
not attempt to validate them** — the data to do that does not exist yet
(see "Benchmark data" in `tasks/todo.md`).

---

## 2. The decision matrix to implement

Per [spec §2.1.3](../spec/02-stages.md#213-decision-matrix). Four branches,
all exiting 0 — a coverage decision is data, not a pipeline error
([principle 8](../CONSTITUTION.md)).

| Condition | `status` | Gated FASTQ | Assembly |
|---|---|---|---|
| `est < HARD_MIN` | `fail` | empty (channel typing) | skipped |
| `HARD_MIN <= est < WARN` | `low_coverage` | full passthrough | **runs** |
| `WARN <= est <= MAX` | `ok` | full passthrough | runs |
| `est > MAX` | `ok` | subsampled | runs |

Boundaries are inclusive at the bottom of each band: exactly `HARD_MIN` is
`low_coverage` (not `fail`), exactly `WARN` is `ok` (not `low_coverage`),
exactly `MAX` is passthrough (not subsampled — existing behaviour, keep).

**The easiest bug to write here** is emitting an empty gated FASTQ on the
`low_coverage` branch by reusing the `fail` branch's file handling. That
would route a sample to the assembler with no reads and produce a confusing
downstream failure instead of a degraded assembly. The `low_coverage`
branch must pass reads through exactly as `ok` does.

---

## 3. Config and parameter migration

In `nextflow.config`, `params.coverage_limits` gains a key and loses one,
per [spec §2.1.1](../spec/02-stages.md#211-per-target-limits):

| Target | `nominal_size` | `hard_min_cov` | `warn_cov` | `max_cov` |
|---|---|---|---|---|
| `animal_mt` | 17000 | 10 | 30 | 300 |
| `plant_pt` | 150000 | 10 | 30 | 500 |
| `plant_mt` | 400000 | 10 | 30 | 300 |

`min_cov` is **retired, not renamed to one of these**. Both new keys carry
part of its old meaning, so a config still setting `min_cov` is ambiguous
and must not be silently absorbed. There is no `nextflow_schema.json` in
this project (only `assets/samplesheet.schema.json`), so schema validation
will not catch it — add an explicit guard where params are materialised at
the top of `main.nf`:

```
# pseudocode, main.nf init
for each target, limits in params.coverage_limits:
    if 'min_cov' in limits:
        error "params.coverage_limits.${target}.min_cov was retired in
               task 35. Set hard_min_cov (terminal floor) and warn_cov
               (advisory floor) instead — see spec §2.1.1."
    if limits.hard_min_cov > limits.warn_cov:
        error "hard_min_cov must not exceed warn_cov for ${target}"
```

The second check matters because the two floors are independently
overridable on the command line (`--coverage_limits.animal_mt.warn_cov 20`),
and an inverted pair would silently collapse the middle band to nothing.

Update the block comment above `coverage_limits` in `nextflow.config`: it
currently describes MIN/MAX as "prior-art defaults, never validated" and
should now state the different provenance of each floor per
[spec §2.1.6](../spec/02-stages.md#216-provenance-of-the-two-floors) — one
observation behind the warn floor, none behind the hard floor.

---

## 4. `bin/coverage_gate.py` (C2)

The sibling-panel split from task 25 is unchanged — it decides *which
bases* feed the estimate, and this task only changes *what the estimate is
compared against*. Do not disturb `map_to_panel` / `merge_intervals` or the
`coverage_basis` fallback logic.

Changes:

1. **CLI:** replace `--min-cov` with `--hard-min-cov` and `--warn-cov`.
   Keep `--max-cov`.
2. **Branching:** implement the four-branch matrix of §2. Structure it so
   the passthrough file-write is shared between the `low_coverage` and `ok`
   passthrough branches rather than duplicated — the §2 bug is much harder
   to write if there is only one copy of that line.
3. **`sample_status.json`:** `min_required` is replaced by
   `hard_min_required` and `warn_threshold`. Every branch emits both, so a
   reader can always see how far the sample sat from each floor. Add
   `below_warn_floor` (boolean) as an explicit convenience flag for report
   templates — deriving it by comparing floats in Jinja invites drift.
4. **Module docstring:** it currently documents the three old outcomes and
   names `status="low_coverage"` as the sub-MIN case. Rewrite it — this
   docstring is the first thing the next reader will trust, and it is
   currently the most misleading text in the file.

`coverage.json` needs no schema change.

---

## 5. Workflow wiring

**`modules/local/coverage_gate.nf`** — pass the two new flags through from
`params.coverage_limits[meta.assembly_target]`. The `stub:` block emits
`{"status": "ok", ...}` and can stay as-is (stub runs must keep exercising
the assembling path), but add a comment noting the vocabulary so the next
person editing it does not invent a fourth status.

**`main.nf`** — the branch condition is currently:

```
ok:     status_json.text.contains('"status": "ok"')
failed: true
```

That substring test excludes `low_coverage`, so degraded samples would be
routed to `COLLATE` and this task's entire intent inverted. The three
statuses share no common prefix, so per
[spec §2.1.3](../spec/02-stages.md#213-decision-matrix) replace it with a
parsed **explicit allowlist**:

```
# pseudocode
status = parseJson(status_json).status
ok:     status in ['ok', 'low_coverage']
failed: true          # 'fail', plus any future status — fail closed
```

Do not substitute a looser substring or prefix match. A status added later
would land on whichever side of the branch its spelling happened to match,
which is the silent mis-routing [rule 18](../CONSTITUTION.md) exists to
prevent. Failing closed (unknown status → no assembly) is the correct
default here: skipping assembly is recoverable, shipping an unvalidated
bundle to Taxodactyl is not.

Rename the branch labels/comments away from "soft-fail" where they now
describe the `fail` status specifically.

---

## 6. Consumer migration — the `low_coverage` redefinition

Every site below reads the old vocabulary. Each must be visited
deliberately; **none of them will fail loudly if missed**.

| Site | Currently | Required |
|---|---|---|
| `bin/coverage_gate.py` docstring + status writes | `low_coverage` = sub-MIN reject | §4 |
| `main.nf` gate branch | substring `"status": "ok"` | §5 allowlist |
| `modules/local/collate.nf` header comment + script comment | "Handles both ok and soft-failed (low_coverage) samples" / "emit full or minimal (low-coverage) bundle" | COLLATE is still a P4 stub, so this is comment-only **plus** a carry-forward note (§6.1) so the real implementation dispatches on three statuses, not two |
| `tests/integration/expected/plant_mt/coverage_bounds.json` | `expected_status: "low_coverage"` meaning *rejected* | §8 |
| `scripts/tests/test_coverage_gate.py` | several tests assert `low_coverage` meaning *rejected* | §7 |
| `tasks/todo.md`, "Integration fixtures" entry | states the fixture soft-fails and cannot support downstream assertions | Update — no longer true after this task; point at task 36 |

### 6.1 Carry-forward for `COLLATE` (P4)

Add to `tasks/todo.md` under the existing `COLLATE` (P4) heading:

> `COLLATE` must dispatch on **three** gate statuses, not two. `fail` gets
> the minimal bundle (no assembly exists). `low_coverage` gets the **full**
> bundle — assembly, annotation and barcodes all exist — with the warning
> carried into `metadata.json` so the per-sample report can surface it and
> the Taxodactyl handoff records that the result is partial
> (spec §2.1.3, §6a.4). Treating `low_coverage` as a failure path would
> discard a real result.

---

## 7. Unit tests (`scripts/tests/test_coverage_gate.py`)

C2 targets **100 % branch coverage** ([rule 14](../CONSTITUTION.md)), and
this task adds a branch. Run via `scripts/pytest.sh`.

**Existing tests that change meaning — check each individually.** This is
the most important part of the test work, because several will keep
passing while asserting something new:

- `test_low_coverage` — 100 kb input ≈ 6× against a 30× floor. Now below
  the 10× hard floor, so it should assert `fail`. Rename to
  `test_below_hard_floor_fails`.
- `test_empty_input` (0 reads) and `test_nominal_size_zero` — both give
  `est = 0`, now `fail`.
- `test_exactly_min_cov_is_passthrough` — rewrite as two tests, one per
  floor boundary (§7 new cases below).
- `test_main_gates_on_target_bases` — its docstring says "Only 4 Mb is
  mitochondrial → 10× → low_coverage" and asserts `status == "low_coverage"`.
  Under the new floors 10× is *exactly* the hard floor, so the assertion
  **still passes** but now means the sample assembles rather than being
  rejected. Left alone this test silently stops testing what it was written
  to test. Move the fixture clearly into one band (e.g. adjust so it lands
  mid-band) and assert the intended outcome explicitly.
- `test_main_subsample_scales_target_estimate` — unaffected, but confirm.

**New cases:**

- Between floors → `status == "low_coverage"`.
- **Between floors emits a non-empty gated FASTQ byte-identical to the
  input.** This is the §2 bug; assert it directly, not via status.
- Below hard floor → empty gated FASTQ (existing behaviour, now on `fail`).
- Boundary: `est == hard_min_cov` exactly → `low_coverage`.
- Boundary: `est == warn_cov` exactly → `ok`.
- `hard_min_required`, `warn_threshold` and `below_warn_floor` present and
  correct in `sample_status.json` on **all four** branches.
- Sibling-split interaction: a sample whose *total* recruited bases clear
  the warn floor but whose *target-assigned* bases land between floors
  reports `low_coverage` — i.e. the split still drives the decision. This
  is the task 25 regression guard and matters more than any other new case,
  because it is where the two features interact.

Follow the boundary-mocking convention of task 27 — fake the
seqkit/seqtk/minimap2 boundary, do not invoke the tools.

---

## 8. Integration reconciliation (minimal here)

`INT-PLANT-01-mt` sits at **28.48×** after task 28. Under the new floors it
clears the 10× hard floor and becomes `low_coverage` — **it now assembles
for the first time**, and stages 7–16 run for `plant_mt` in CI where they
previously did not.

Two consequences, and only the first belongs to this task:

1. **CI budget.** METAFLYE onward now runs for a third sample, against the
   ≤ 60-minute wall-clock budget in
   [spec §5](../spec/05-test-data.md). Measure the new total runtime of
   `-profile integration` and record it in this task's completion notes.
   If it breaches the budget, say so — do not silently trim assertions to
   fit; that is task 36's problem to solve with the measurement in hand.

2. **Downstream assertions** for the newly-reachable `plant_mt` path —
   **task 36**, not this one.

In this task, update `tests/integration/expected/plant_mt/coverage_bounds.json`:

- Keep `expected_status: "low_coverage"` — but rewrite `_comment`, which
  currently explains the value as a soft-fail against a 30× MIN. As written
  it would leave a future reader with exactly the wrong model, while the
  assertion goes on passing.
- Bounds `min/max_estimated_cov` (15–30) still bracket 28.48× and can stay.

Add one assertion to `tests/integration/assertions.sh` that pins the *new*
semantics rather than the string: for a `low_coverage` sample, assert
`gated.fastq.gz` is **non-empty**. Under the old meaning it was empty by
construction, so this single check distinguishes the two regimes and would
have caught the §2 bug. Leave `ASSEMBLING_SAMPLES` alone — task 36 owns
that.

---

## 9. Acceptance criteria

- `nextflow.config` carries `hard_min_cov` / `warn_cov` per target; a
  config setting `min_cov`, or one with `hard_min_cov > warn_cov`, aborts
  the run with a message naming the spec section.
- C2 implements the four-branch matrix; boundaries resolve as specified in
  §2; the between-floors branch passes reads through unmodified.
- `main.nf` branches on a parsed explicit allowlist and fails closed on an
  unrecognised status.
- `pytest` passes with C2 at 100 % branch coverage; every test named in §7
  has been individually re-examined, not just left green.
- `flake8` clean on `bin/coverage_gate.py` (`claude` venv, per CLAUDE.md).
- `nextflow run . -profile stub -stub-run` passes.
- `-profile integration` passes, with `INT-PLANT-01-mt` reaching METAFLYE;
  total wall-clock recorded against the 60-minute budget.
- Every row of §6's migration table is closed, including the `todo.md`
  edits.

---

## 10. Out of scope

- **Validating either floor.** Both are provisional placeholders
  (spec §9 item 2), and the benchmark data to sweep them does not exist.
  Do not tune them against the integration fixtures — that is circular
  ([spec §5.1](../spec/05-test-data.md#51-the-integration-fixtures-are-correctness-tests-not-benchmarks)).
- **`plant_mt` downstream integration assertions** — task 36.
- **Report rendering** of the warning. The Jinja report work is P4a's
  `COLLATE`/`RUN_REPORT` boilerplate; this task only guarantees the fields
  exist in `sample_status.json` and adds the §6.1 carry-forward.
- **Replacing `INT-PLANT-01-mt`** with a deeper fixture. Still wanted
  (`tasks/todo.md`) — a sample that only clears the *hard* floor cannot
  validate the `plant_mt` arm, it can only exercise it.
