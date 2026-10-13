# Task 36 — `plant_mt` integration coverage now that the fixture reaches assembly

**Phase:** P4a (integration-test completion), gating P4b.
**Goal:** Turn the `plant_mt` integration path from "runs but is unasserted"
into real end-to-end coverage, now that `INT-PLANT-01-mt` clears the hard
coverage floor and reaches METAFLYE for the first time.

**Status:** blocked on task 35. Do not start before it lands — this task's
entire premise is the behaviour change task 35 introduces.

**Related:** consumes the runtime measurement from task 35 §8. Closes the
"Integration fixtures" entry in `tasks/todo.md` *partially* — the deeper
fixture that entry asks for is still wanted, for different reasons (§5).

---

## 1. Overview

Integration testing in this project runs real tools against small real ONT
fixtures, once nightly, and asserts on the outputs
([spec §5](../spec/05-test-data.md)). There are three fixtures, one per
assembly target.

One of them, `INT-PLANT-01-mt`, has never been tested past stage 6. Its
recruited read pool is about 70 % chloroplast carry-over, so once task 25
taught the coverage gate to estimate from mitochondrial bases only, the
sample landed at 28.48× against a 30× floor and was rejected before
assembly. That was the *correct* verdict, but it left the entire `plant_mt`
branch — assembly, contig binning, the multi-contig `emit: all` path,
sibling-organelle discrimination, barcode extraction — covered by unit
tests alone, with no real-tool exercise at all.

Task 35 splits that floor in two. At 28.48× the sample now sits between the
hard floor (10×) and the warn floor (30×), so it assembles and flows
through every remaining stage. The pipeline stages are already running in
CI; what is missing is anything checking that their output is sensible.

This task writes those checks. The framing matters, and §2 is the whole
difficulty: **this fixture is a degraded sample, and the pipeline itself
says so.** Its assertions must pin the behaviour we actually want from a
degraded input — that the pipeline completes, stays self-consistent, and
labels the result honestly — without quietly enshrining a partial assembly
as the standard the `plant_mt` arm is held to.

---

## 2. The framing constraint (read before writing any assertion)

`INT-PLANT-01-mt` is flagged `low_coverage` by the gate. Every downstream
number it produces is a degraded number.

That makes two kinds of assertion legitimate and one kind actively harmful:

**Legitimate — structural and self-consistency checks.** Does the stage
emit its expected files? Is the GFF internally consistent with the FASTA?
Does `bin_metadata.json` classify contigs into the schema's categories?
Does the sibling-organelle fraction stay in a sane range? These test that
the code is correct, and they hold regardless of depth.

**Legitimate — floor assertions, framed as floors.** "At least one
mitochondrial contig is selected", "at least one barcode locus is
recovered". Written as minima with generous headroom, these catch a real
regression (the arm silently producing nothing) without asserting that the
current degraded output is *right*.

**Harmful — tight bounds derived from this run's output.** Recording the
assembly length, contig count or locus count this fixture happens to
produce and asserting near-equality. That is
[spec §5.1](../spec/05-test-data.md#51-the-integration-fixtures-are-correctness-tests-not-benchmarks)'s
circularity trap with an extra twist: the numbers would come from a sample
the gate has already labelled as not trustworthy for completeness. A future
genuine improvement to recruitment or assembly would then fail CI for
producing a *better* assembly.

`tests/integration/expected/plant_mt/{assembly,bin}_bounds.json` already
exist from earlier work and are currently retained-but-unasserted. **Treat
their contents as suspect, not as a starting point** — they were written
before task 25's sibling split and task 28's panel change, against a
different recruited pool. Re-derive every bound from the current run and
widen it deliberately; do not simply switch them on.

---

## 3. Scope

### 3.1 Re-enable the sample

Add `INT-PLANT-01-mt` to `ASSEMBLING_SAMPLES` in
`tests/integration/assertions.sh`. That single change activates the
existing per-stage assertion blocks (METAFLYE, BANDAGE_NG, BIN_TARGET,
BLAST_VALIDATE, EXTRACT_BARCODES, ANNOTATE) for this sample. Work through
the resulting failures one stage at a time rather than pre-emptively
loosening bounds.

Note `ANNOTATE` on `plant_mt` is CDS-only (`non_cds_tool: 'none'`), so its
summary carries `status: "ok_cds_only"` and zero tRNA/rRNA counts — expected,
not a failure ([spec §8 item 3](../spec/07-open-questions.md#8-remaining-open-questions)).

### 3.2 Bounds files

Re-derive `assembly_bounds.json` and `bin_bounds.json` per §2. Each bound
carries a `_comment` recording that it was set from a `low_coverage`
fixture and is a floor, not a target.

### 3.3 The degraded-status assertion

Assert that the degraded state is *visible*, not just tolerated:
`sample_status.json` reports `low_coverage`, and whatever downstream
artefacts carry the flag (per task 35 §4's `below_warn_floor`) agree with
it. A pipeline that assembles a degraded sample and then presents it as
clean is the exact failure [principle 7](../CONSTITUTION.md) forbids, and
this fixture is the only place we can test that end-to-end.

### 3.4 CI budget

Task 35 §8 records the new `-profile integration` wall-clock. If it breaches
the 60-minute budget in [spec §5](../spec/05-test-data.md), resolve it here.
Options, in preference order:

1. Accept it if there is headroom — simplest.
2. Trim the fixture's read count. **Risky and probably self-defeating:**
   the sample is at 28.48×, and the hard floor is 10×; subsampling it
   further walks it toward `fail` and back to testing nothing. Whatever is
   removed must leave a comfortable margin above the hard floor, and the
   margin must be stated in the fixture's `_comment`.
3. Split the nightly job per target so each fixture runs in its own
   parallel CI job.

Do not resolve a budget breach by removing assertions.

---

## 4. Acceptance criteria

- `INT-PLANT-01-mt` is in `ASSEMBLING_SAMPLES`, and `-profile integration`
  passes end-to-end with every stage from METAFLYE to ANNOTATE asserted.
- Bounds are floors with headroom, each carrying a `_comment` recording its
  degraded provenance; none is a near-equality on this run's output.
- The `low_coverage` flag is asserted as visible downstream, not merely
  present at the gate.
- Total integration wall-clock is inside the 60-minute budget, or the
  breach is resolved by §3.4 and the resolution recorded.
- `tasks/todo.md`'s "Integration fixtures" entry is rewritten: the
  `plant_mt` branch is no longer uncovered, but a deeper fixture is still
  wanted per §5.

---

## 5. What this task does *not* achieve

It gives the `plant_mt` arm **correctness** coverage, not **validation**.

A fixture that only clears the hard floor can show the arm runs end-to-end
and stays self-consistent. It cannot show the arm produces a good
mitogenome, because the input is not good enough for that question to be
meaningful. Two things therefore remain open regardless of how this task
goes:

- The deeper `plant_mt` fixture in `tasks/todo.md` — a sample clearing the
  **warn** floor, needed before anything about `plant_mt` assembly quality
  can be claimed.
- P4b's annotation benchmark
  ([spec §9 item 12](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)),
  which compares annotators against a curated ground truth and needs a
  sample whose assembly is trustworthy. Per
  [spec §6](../spec/06-phases.md), this task relieves P4b's fixture gate —
  the arm can now be *built* against a running fixture — but does not close
  it.

State this plainly in the completion notes. The risk this task creates is a
false sense of `plant_mt` coverage: CI turning green on a sample the
pipeline itself flags as degraded is easy to misread as the arm being
finished.
