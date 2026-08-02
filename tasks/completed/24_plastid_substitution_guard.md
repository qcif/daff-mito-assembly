# Task 24 — Guard C4's `target.fasta` substitution behind C3 selection

**Phase:** P3 (from [spec §6](../../spec/06-phases.md)).
**Goal:** Close a silent false-positive path in
[`bin/bin_target.py`](../../bin/bin_target.py): on the `plant_pt` canonical
branch, C4's canonicalised `path1.fasta` **unconditionally** overwrites
`target.fasta`, regardless of whether C3 selected any plastid contig. A
sample whose assembly contains nothing recognisable as a plastid can
therefore still emit a confident ~150 kb plastome.

**Prerequisite:** [task 23](23_bin_target_recalibration.md) — it touches
the same selection path, and its fixes are what make C3's `plant_pt`
selection meaningful in the first place. Sequence 22 before 23.

**Spec basis:**
[§3.6 step 5](../../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation),
amended to state the precondition this task implements.

## 1. The defect

[`bin/bin_target.py`](../../bin/bin_target.py) on the `plant_pt` branch:

```
if result.branch == 'canonical':
    shutil.copyfile(iso_dir / 'path1.fasta', args.out_target)
```

There is no reference to `primaries` in that condition. C4 keys entirely
off the **graph topology** (a 3-edge `assembly_graph.gfa`); C3's
contig-level evidence — reference homology, coverage, ORF — is discarded
whenever the graph happens to have three edges.

This is currently masked because it produces the *right* answer on the
integration fixture for the *wrong* reason. `INT-PLANT-01-pt` reports:

```json
{ "n_target_selected": 0,
  "contigs_selected": [],
  "plastid_canonicalisation": { "branch": "canonical",
                                "substitution_applied": true } }
```

C3 selected **nothing**, yet `target.fasta` is a correct 155 kb
plastome, and
[`tests/integration/assertions.sh`](../../tests/integration/assertions.sh)
passes it — the assertion only checks `target.fasta` is non-empty and
equals `path1.fasta`. The `plant_pt` binning path is, in effect,
currently untested.

## 2. Why it matters beyond the test

Three edges in a Flye graph is a *structural* observation, not a
taxonomic one. A 3-edge graph arises from any assembly with one
high-depth repeat between two single-copy segments. The branch does not
verify that the edges are plastid — that is precisely C3's job, and its
verdict is being thrown away.

Consequences, in order of seriousness:

1. **A misdeclared sample gets a confident answer.** An `animal_mt`
   sample submitted as `plant_pt`, or a plant sample with no usable
   plastid recovery, yields a 3-edge graph and ships a ~150 kb
   `target.fasta` that BLAST_VALIDATE, ANNOTATE, and MINIPROT_EXTRACT
   all consume as the plastome. Directly contrary to
   [principle 7](../../CONSTITUTION.md) (a degraded sample must never
   produce output indistinguishable from a confident negative) and
   [rule 18](../../CONSTITUTION.md) (a subtly wrong answer that ships as a
   "result" is worse than a loud failure).
2. **The `no_assembly` negative state is unreachable on `plant_pt`.**
   [Principle 7](../../CONSTITUTION.md) requires three distinguishable
   negative signals. On the canonical branch `target.fasta` is never
   empty, so `plant_pt` can never report `no_assembly` however poor the
   recovery.
3. **C3's `plant_pt` path has no test coverage in integration.** Any
   regression in `plant_pt` binning is invisible while the substitution
   masks it. Task 23's `n_target_selected >= 1` assertion is what
   exposes it.

## 3. Change

Make the substitution conditional, and record the disagreement:

```
# pseudocode
canonical = result.branch == 'canonical'
selected  = len(primaries) > 0

if canonical and selected:
    substitute path1.fasta -> target.fasta      # as today
elif canonical and not selected:
    keep target.fasta empty (C3's verdict stands)
    keep plastid_isoforms/ as diagnostics
    metadata: substitution_applied = false
              substitution_withheld_reason = 'no_c3_selection'
else:
    unchanged from current behaviour
```

The isoform files stay on disk in the withheld case — they are useful
diagnostics for an operator working out *why* the sample failed, and
[principle 7](../../CONSTITUTION.md) favours more signal, not less. Only
the substitution into the contractual output is withheld.

### 3.1 Report surface

The withheld case is a QC signal an operator must see, not a silent
metadata field. Add to the per-sample report's "Assembly quality
assessment" section
([§6a.2](../../spec/06a-reports.md)), alongside the existing edge-count
signal: *"plastid graph is canonical (3 edges) but no contig passed
target binning — assembly structure looks plastid-like, taxonomic
evidence does not support it. Manual review recommended."*

Also surface the converse as provenance in the normal case: that
`target.fasta` is C4's `path1`, not the raw C3-selected contig. An
auditor reading `contigs_selected: ["contig_2"]` today has no indication
the emitted sequence is a reconstructed path rather than that contig.

### 3.2 Interaction with `COLLATE` (P4)

`COLLATE` dispatches on the coverage-gate marker and, in P4, on the
`no_assembly` state. An empty `target.fasta` from the withheld branch
must route to the `no_assembly` bundle path, not to the full bundle with
empty downstream outputs. This task does not implement `COLLATE` — flag
it in the `COLLATE` task's inputs so the state is handled when that
lands.

## 4. Tests

**Unit** — `scripts/tests/test_plastid_canonicalise.py` /
`test_bin_target.py`, per [rule 14](../../CONSTITUTION.md):

| # | Scenario | Expected |
|---|---|---|
| 1 | Canonical 3-edge graph, C3 selected ≥ 1 contig | Substitution applied; `target.fasta == path1.fasta`; `substitution_applied true` |
| 2 | Canonical 3-edge graph, C3 selected 0 contigs | `target.fasta` empty; `plastid_isoforms/` present; `substitution_applied false`; `substitution_withheld_reason` set |
| 3 | `resolved_circle` branch, C3 selected 0 | Unchanged from current behaviour; no isoform dir |
| 4 | `non_canonical` branch, C3 selected ≥ 1 | C3's own selection is emitted; no substitution |

Case 2 is the regression test for this defect.

**Integration** — covered by task 23's `n_target_selected >= 1`
assertion on all three fixtures. Once task 23 lands, `INT-PLANT-01-pt`
selects contig_2 on C3's own evidence, so the substitution proceeds via
the guarded path and the existing
`target.fasta == plastid_isoforms/path1.fasta` assertion still holds —
but now for the right reason. Add an assertion that
`substitution_applied == true` **and** `n_target_selected >= 1`
together, so the two can never again diverge unnoticed.

## 5. Exit criteria

- `bin_metadata.json` on `plant_pt` never reports
  `substitution_applied: true` with `n_target_selected: 0`.
- Withheld case emits empty `target.fasta`, retains `plastid_isoforms/`,
  records the reason, and exits 0 (a data outcome, not a Nextflow error,
  per [principle 8](../../CONSTITUTION.md)).
- Report renders the "canonical structure, no taxonomic support" warning.
- `INT-PLANT-01-pt` integration passes with both
  `n_target_selected >= 1` and `substitution_applied == true`.
- 100 % branch coverage on the substitution decision;
  `scripts/pytest.sh` green; flake8 clean.
- [spec §3.6 step 5](../../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation)
  already carries the precondition — verify the implementation matches it.

## Outcomes

Completed 2026-08-03. Implemented as specified in §3; no deviation from
the brief.

### What changed

- [`bin/bin_target.py`](../../bin/bin_target.py) — the substitution
  condition gained a second conjunct. It was C4's graph verdict alone
  (`result.branch == 'canonical'`); it is now
  `canonical and bool(primaries)`, where `primaries` is C3's list of
  selected target contigs, so a non-empty list means C3 found
  reference-homology evidence for a plastid. The reconstructed `path1`
  can overwrite `target.fasta` only when the graph *looks* plastid-like
  **and** the taxonomic evidence agrees it is one. Three new constants
  (`TARGET_SOURCE_C3`, `TARGET_SOURCE_C4_PATH1`,
  `SUBSTITUTION_WITHHELD_NO_SELECTION`) keep the metadata vocabulary out
  of the branch logic. The `shutil.rmtree(iso_dir)` cleanup is now
  guarded by `not canonical`, so the withheld case keeps its isoforms as
  operator diagnostics while every other branch still clears upstream
  residue.
- **New `target_source` metadata field** (§3.1's "converse as
  provenance"): `c3_selected_contigs` by default on every target,
  `c4_plastid_path1` when the substitution applied. An auditor reading
  `contigs_selected: ["contig_2"]` can now tell whether `target.fasta`
  is that contig or C4's reconstructed path. Recorded in
  [spec §3.6 step 5](../../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation).
- [`tests/integration/assertions.sh`](../../tests/integration/assertions.sh)
  — the canonical branch now asserts `substitution_applied == true`
  **and** `n_target_selected >= 1` together, plus
  `target_source == c4_plastid_path1`.
- [`spec/03-organelles.md`](../../spec/03-organelles.md) §3.6 — added the
  concrete metadata field names to step 5 and a "withheld-substitution
  QC signal" paragraph alongside the existing edge-count signal.
- [`spec/06a-reports.md`](../../spec/06a-reports.md) §6a.2 — the
  "Assembly quality assessment" row now names the three plastid QC
  signals it must carry.

### Deferred (not implementable in P3)

§3.1's report rendering and §3.2's `COLLATE` routing both land in P4 —
neither the per-sample report nor `COLLATE`'s `no_assembly` dispatch
exists yet. Both are specified above (spec §3.6 / §6a.2) and carried as
entries in [`tasks/todo.md`](../todo.md) so the drafting of those tasks
inherits them.

### Verification

- `scripts/pytest.sh` — 95 passed. `bin/bin_target.py` at 99 %
  branch coverage; the only uncovered line is the `if __name__ ==
  '__main__'` guard, so the substitution decision itself is at 100 %.
  The 4 `test_coverage_gate.py` failures are the pre-existing missing-
  `seqkit` issue owned by [task 27](../27_unit_test_boundary_mocking.md).
- Five unit tests cover the §4 matrix — cases 1–4 plus a fifth for the
  stale-isoform cleanup branch that the guard's restructuring exposed.
  Case 2 (`test_plant_pt_canonical_without_selection_withholds`) is the
  regression test; it fails against the pre-task code.
- `-profile integration` (2026-08-03, full re-run, 47 processes,
  25 m 30 s). `bash tests/integration/assertions.sh` — all assertions
  passed, including the two new ones:
  `INT-PLANT-01-pt substitution_applied with n_target_selected=1` and
  `target_source=c4_plastid_path1`.
- `flake8` clean on both changed Python files.

### Note on the brief

§Prerequisite reads "Sequence 22 before 23" — a typo for "sequence 23
before 24". Task 23 was already complete, so the intended ordering held.
