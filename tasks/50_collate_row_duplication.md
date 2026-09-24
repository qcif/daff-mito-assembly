# Task 50 — `ch_ok_inputs` emits duplicate rows per sample before `COLLATE`

## Overview

`COLLATE` (C6, task 42) is supposed to receive exactly one input tuple
per sample — one bundle in, one bundle out (CONSTITUTION hard
constraint 3: one sample row → one organelle assembly). During task
44's first real `-profile integration` run against the pinned
production containers, `main.nf`'s `ch_ok_inputs` channel — the 13-way
join chain at [main.nf:251-285](../main.nf#L251-L285) that assembles
`COLLATE`'s input tuple from every upstream stage's output — was
observed emitting **more than one row for every assembling sample**:
4 rows for `INT-ANIMAL-01`, 4 for `INT-PLANT-01-mt`, and 9 for
`INT-PLANT-01-pt`.

This blocks every sample's `report.html` from rendering, and therefore
blocks `RUN_REPORT` and the whole tail of the pipeline.

**Task 44 did not introduce this.** `ch_ok_inputs` and every join
feeding it predate task 44 untouched — task 44 only added one new
`.join(BIN_TARGET.out.metadata, by: 0)` earlier in the workflow, at the
`ORGANELLE_MAP` call site, which is 1:1 by construction (`BIN_TARGET`
emits `bin_metadata.json` exactly once per sample). This task exists
because task 44's acceptance criteria needed a live integration run to
validate the new organelle map renderer, and that run is what surfaced
this pre-existing defect.

It has apparently never been caught before: task 47 §5.4 records that
`tests/integration/output/` is stale (pre-task-43a) and `report.html`
is zero bytes for all three samples in the committed fixture snapshot,
meaning no one has watched a full `-profile integration` run reach
`COLLATE`/`RUN_REPORT` with real container output since early in the
project.

## What's confirmed

- Instrumentation used: a `.view { "DEBUG_ROW size=${it.size()} : ${it}" }`
  inserted directly after `.join(BIN_TARGET.out.isoforms, by: 0,
  remainder: true)` (the last join in the chain) and before the row's
  final `.map` in `main.nf`. **This instrumentation has since been
  removed** — step 1 below re-establishes it.
- Every emitted row has the correct flat width (21 elements, matching
  the closure's declared arity) — this is **not** an arity mismatch.
  It is the same correct row shape, repeated.
- Row counts per sample: `INT-ANIMAL-01` ×4, `INT-PLANT-01-mt` ×4,
  `INT-PLANT-01-pt` ×9.
- Downstream, `COLLATE` was observed to abort the whole Nextflow
  session with a Groovy `MissingMethodException` on a different closure
  once real (non-cached) data reached it. This is plausibly a
  consequence of the duplication rather than a separate defect —
  re-check it once duplication is fixed; it may simply disappear.

## Why the obvious hypothesis is wrong

The counts are exact squares (2², 2², 3²), which *looks* like a
relational cross-product. **Do not start from that assumption** —
Nextflow's `join` operator cannot produce one.

`join` buffers items by key and, for each key, emits one combined tuple
per *matched set*, popping one item from each input channel. With
duplicate keys the emitted count for a key is therefore `min()` of the
per-channel counts, never their product. Two channels each emitting a
sample 2× yield 2 joined rows, not 4.

Two consequences worth holding onto while investigating:

- Under standard `join` semantics, 4 rows out of the chain requires
  **every one of the 13 input channels** to have emitted that sample
  at least 4 times. That would mean the fan-out originates at or above
  `ch_gated.ok` and propagates through the entire assembly chain — in
  which case `METAFLYE`, `BIN_TARGET` etc. each ran 4× for that sample
  and the duplication is visible in the process invocation counts.
- If the invocation counts show each process running **once** per
  sample, then the duplication is being introduced by an operator, not
  a process, and the `min()` reasoning above is being defeated by
  something specific — most likely `remainder: true` on the final join
  (which changes emission semantics and is the one non-standard join in
  the chain), or a value-vs-queue channel being replayed. Either way
  that observation sharply narrows the search.

The square-number pattern may simply be coincidence. Let the
measurement decide.

## Investigation plan

1. **Re-establish the measurement.** Re-add the `.view()` after the
   final join in `ch_ok_inputs` and confirm the 4/4/9 counts still
   reproduce on a clean `-profile integration` run (delete `./work/`
   first — do not `-resume` into stale task data). Reducing this to a
   per-key tally is easier to read than raw rows, e.g.:

   ```groovy
   ch_ok_inputs.map { it[0].sample_id }
       .groupTuple( by: 0 )   // or simply .count() per branch
       .view { "ROWCOUNT ${it}" }
   ```

2. **Localise with process invocation counts — do this before any
   bisection.** Run with `-with-trace` (no trace config exists in
   `nextflow.config`/`conf/`, so pass it on the command line) and count
   executions per process per sample:

   ```bash
   nextflow run . -profile integration -with-trace trace.txt
   awk -F'\t' 'NR>1 {print $4}' trace.txt | sort | uniq -c | sort -rn
   ```

   Any process showing >1 invocation per assembling sample localises
   the fan-out to the *first* such stage in the DAG. This single step
   likely replaces most of the 13-way bisection below.

3. **Bisect only if step 2 is inconclusive.** Add a per-key row-count
   `.view()` after each `.join(..., by: 0)` in turn and find where the
   count first exceeds 1: `NANOPLOT_RAW.out.reports`,
   `NANOPLOT_CLEAN.out.reports`, `RECRUIT.out.stats`,
   `BIN_TARGET.out.binned`, `BLAST_VALIDATE.out.validated`,
   `EXTRACT_BARCODES.out.barcodes`, `ANNOTATION_SCORING.out.annotation`,
   `ORGANELLE_MAP.out.map`, `ch_graph_png`, `BIN_TARGET.out.metadata`,
   `ch_assembly_info`, `ch_genetic_code`, `BIN_TARGET.out.isoforms`.

   Prior suspicions worth checking early, given `plant_pt`'s higher
   count: any stage emitting per-contig rather than per-sample tuples
   for a multi-contig target — `plant_mt`'s "up to 20 contigs"
   `emit: all` path (task 23) and `plant_pt`'s primary+secondary pair
   (task 24). Also note that `ch_cds` is consumed four times
   ([main.nf:202-244](../main.nf#L202-L244): `EXTRACT_BARCODES`,
   `ANNOTATE`, `ANNOTATION_SCORING`, `ch_genetic_code`) and
   `BLAST_VALIDATE.out.validated` and `BIN_TARGET.out.metadata` twice
   each — confirm each fork behaves as expected.

4. **Establish root cause.** Determine whether the offending stage's
   process declares one output tuple per sample (correct) or one per
   contig/file (incorrect for a sample-keyed join), or whether the
   process is correct and the **operator** is the defect (a value vs
   queue channel being replayed, a `remainder: true` interaction, or a
   channel that should have been `.groupTuple()`-reduced to one row per
   sample before joining).

5. **Fix at whichever layer is actually wrong** — a process/subworkflow
   emit that needs reducing to one row per sample, or a join that needs
   to key on more than `sample_id`.

6. **Re-check the `COLLATE` `MissingMethodException`.** Confirm it is
   resolved by the fix; if it survives, scope it as a further finding
   in this task's Outcomes.

## Constraints

- Do not change what any stage's output *means* — only how many tuples
  it contributes to the sample-keyed join. A fix must leave every
  currently-correct multi-file output (e.g. `plant_mt`'s multi-contig
  `target.fasta`/`secondaries.tsv`, which is intentionally
  multi-record *within* one file, not one row per contig) exactly as it
  is.
- CONSTITUTION hard constraint 3 (one sample row → one organelle
  assembly) is the acceptance bar: `COLLATE` must receive exactly one
  tuple per sample, always.
- CONSTITUTION principle 8 (cross-sample failure isolation): whatever
  the fix, one sample's fan-out must not be able to abort the session
  for the others — if the crash in step 6 turns out to be independent,
  it needs its own guard.

## Tests

- A clean `-profile integration` run (real containers, real fixtures,
  `./work/` deleted, no `-resume`) showing exactly one `COLLATE`
  invocation per sample in `trace.txt`, and a non-empty `report.html`
  for all three samples plus `run-report.html`.
- `tests/integration/assertions.sh` already checks `report.html`
  content per sample ([assertions.sh:186](../tests/integration/assertions.sh#L186)
  onwards) but has no assertion on *invocation count* — a silent
  multi-invocation bug republishes the same path and leaves those
  checks green. Add an assertion that would have caught this, e.g.
  parsing `trace.txt` for exactly one `COLLATE` task per entry in
  `ASSEMBLING_SAMPLES`.
- Unit tests only if the fix lands in `bin/*.py` (CONSTITUTION rule 14:
  100% branch coverage). A pure `main.nf` channel fix is covered by the
  integration assertion above; also confirm `-profile stub` still
  passes, since the stub run exercises the same channel topology.

## Acceptance checklist

- [ ] Duplication reproduced and measured on a clean integration run
- [ ] Fan-out localised to a specific stage or operator, with evidence
      recorded in Outcomes
- [ ] Root cause identified (process emit vs channel operator)
- [ ] Fix applied at the correct layer, within the constraints above
- [ ] Exactly one `COLLATE` invocation per sample confirmed via trace
- [ ] `COLLATE` `MissingMethodException` resolved, or scoped as a
      separate finding
- [ ] Integration assertion added guarding invocation count
- [ ] `-profile stub` still green
- [ ] Non-empty `report.html` for all three fixtures rendered, and the
      path given to the user
