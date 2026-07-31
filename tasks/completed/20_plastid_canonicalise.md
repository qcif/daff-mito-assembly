# Task 20 — P3 stage 10 helper: `plastid_canonicalise` implementation (C4)

**Phase:** P3 (from [spec §6](../../spec/06-phases.md)).
**Goal:** Implement the plastid quadripartite canonicalisation
algorithm — specified in
[spec/plastid-canonicalisation.md](../../spec/plastid-canonicalisation.md)
— as `bin/plastid_canonicalise.py` (custom-logic component C4 in
[spec §2.2](../../spec/02-stages.md#22-custom-logic-components)); wire
the `isoforms` emit on `BIN_TARGET`; teach `bin/bin_target.py` to
substitute `path1` for `target.fasta` on the `plant_pt` canonical
branch; record the canonicalisation outcome in `bin_metadata.json`.

C4 is **not a new pipeline stage** — it is a helper library invoked
from within C3 on the `plant_pt` branch of stage 10 (BIN_TARGET). It
shares the `neoformit/daff-wf5-bin-target` container built in
[task 18](18_bin_target.md) and emits additional outputs via the
`emit: isoforms` block on the existing BIN_TARGET module. No workflow
topology change.

---

## Sole algorithm reference

[spec/plastid-canonicalisation.md](../../spec/plastid-canonicalisation.md)
is a self-sufficient specification: prose description, edge
classification rules, path construction rules, output format, CLI arg
surface, library entry-point contract, a test-case matrix, and an
explicit list of implementation choices left open to the implementer.

**Work solely from that document.** It is the only permitted algorithm
reference for this task. Specifically:

- Every detail needed to write `bin/plastid_canonicalise.py` is in the
  spec. If something genuinely is not, **stop and file a follow-up to
  extend the spec** rather than inferring the behaviour or searching
  the repository for anything resembling prior art.
- Do not vendor, import, or adapt any third-party implementation of
  this algorithm, whether found in this repository or elsewhere. The
  deliverable is original code written from the specification.
- The spec's §10 ("implementation-choice boundaries") marks the
  decisions that are deliberately yours (sequence class, regex shape,
  tie-break rule, `Result` type, empty-sequence handling). Any choice
  inside the stated range is acceptable; document the one you make in
  a comment.
- Where this task file and the spec disagree on any algorithmic or
  interface detail, **the spec wins** — this file describes the
  surrounding wiring, not the algorithm.

---

**Prerequisites (both hard-blockers):**

1. [spec/plastid-canonicalisation.md](../../spec/plastid-canonicalisation.md)
   is committed and final. Without it this task cannot start — there
   is no other permitted algorithm reference.
2. [task 18 — BIN_TARGET](18_bin_target.md) is real, so `plant_pt`
   samples flow a real `assembly_graph.gfa` and a real primary
   `target.fasta` through the module. This task additionally tightens
   the plant_pt `bin_bounds.json` that task 18 sets loose (§4 below).

**Exit criteria:**

- For `INT-PLANT-01-pt`:
    - `results/INT-PLANT-01-pt/bin_target/plastid_isoforms/{path1,path2}.fasta`
      both present (3-edge canonical case expected on the SRR-derived
      plastid fixture), each ≥ 100 kb.
    - `target.fasta` **is** `path1.fasta` content (byte-identical
      after substitution); total length ≈ LSC + 2·IR + SSC
      (~145–160 kb).
- `bin_metadata.json.plastid_canonicalisation` records the full
  `Result` payload from
  [spec §4.2](../../spec/plastid-canonicalisation.md): branch
  (`canonical` / `resolved_circle` / `non_canonical`), edge count, the
  assigned LSC / IR / SSC edge names on `canonical`, path lengths, the
  non-canonical reason string where applicable, and whether the
  `target.fasta` substitution was applied.
- `animal_mt` and `plant_mt` samples produce **no**
  `plastid_isoforms/` directory (guard:
  `meta.assembly_target == 'plant_pt'`).
- `tests/integration/expected/plant_pt/bin_bounds.json`
  `target_min_bp` tightened from `70000` → `140000` (see §4).
- `bin/plastid_canonicalise.py` has unit tests covering the full
  test-case matrix in
  [spec §9](../../spec/plastid-canonicalisation.md) using synthetic
  GFAs, with 100 % branch coverage on the classifier and
  canonicalisation functions per
  [CONSTITUTION.md rule 14](../../CONSTITUTION.md).
- Integration assertion added to
  [`tests/integration/assertions.sh`](../../tests/integration/assertions.sh)
  — plant_pt-only: `plastid_isoforms/path1.fasta` exists and matches
  `target.fasta` byte-for-byte.
- Both `-profile stub -stub-run` (fast CI) and `-profile integration`
  (nightly) go green end-to-end.

**Not in scope:**

- Container rebuild — C4 reuses the `neoformit/daff-wf5-bin-target`
  image built in [task 18](18_bin_target.md). `biopython` is already
  installed in it; no additional deps.
- Downstream stage changes (`BLAST_VALIDATE`, `ANNOTATE`,
  `MINIPROT_EXTRACT`) — they continue to consume `target.fasta`
  unchanged. The plastid-specific logic is confined to C3 + C4 per
  [spec §3.6](../../spec/03-organelles.md).
- Rendering both isoforms in `ORGANELLE_MAP` — deferred to the
  ORGANELLE_MAP task drafting, tracked in [tasks/todo.md](../todo.md).
- Handling of GFA `L` (link) lines — canonicalisation reads only `S`
  lines (edge sequence + depth). Link orientation is BandageNG's
  concern.

**Cross-cutting rules (from [spec §1a](../../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- No container change — inherit task 18's SHA pin.
- No host tools; Python + `biopython` come from the C3/C4 shared
  container.
- `bin/plastid_canonicalise.py` is auto-staged onto every process
  `PATH` by Nextflow's `bin/` convention — importable by C3 because
  both scripts sit under `bin/`.

---

## 1. `bin/plastid_canonicalise.py` (new — custom logic C4)

Python + `biopython` (already in the C3 container). The algorithm,
input schema, output format, CLI arg surface and library entry-point
contract are all defined in
[spec/plastid-canonicalisation.md](../../spec/plastid-canonicalisation.md)
— implement to that document, not to this summary.

Required shape of the module (structural requirements only; see the
spec for behaviour):

- **A GFA parsing helper** that reads `S` lines only and yields, per
  edge, its name, sequence, depth and length. Depth comes from the
  `dp:f:` / `DP:f:` tag; a missing tag means depth `0.0`
  (spec §3.2, §10).
- **A classification + canonicalisation entry point** named
  `canonicalise_plastid(gfa_path, outdir=".") -> Result`, matching the
  signature and side-effect contract in
  [spec §8](../../spec/plastid-canonicalisation.md). This is the function
  C3 imports and the function the unit tests exercise.
- **A `Result` type** carrying exactly the fields listed in
  [spec §4.2](../../spec/plastid-canonicalisation.md), including
  `non_canonical_reason`.
- **A thin CLI wrapper** (`main()` + `if __name__ == "__main__"`)
  exposing the positional GFA path plus `--outdir` and `--json-out`
  per [spec §7](../../spec/plastid-canonicalisation.md), printing a
  one-line human-readable summary to stdout.

Design notes:

- **`biopython` only.** No new deps — the C3 container from task 18
  already ships `biopython` for ORF search.
- **Pure library + thin CLI.** `canonicalise_plastid()` is the tested
  interface; the CLI is a wrapper for standalone debugging. C3
  imports the function directly.
- **Module import must be side-effect free** so C3 can import it
  inside a conditional branch (spec §8).
- **Total classifier.** Every branch the spec enumerates — including
  the degenerate LSC/IR-collision and depth-tie cases — must return a
  `Result` with the specified reason string rather than raise. Only
  genuinely corrupt input raises (spec §8).

## 2. `bin/bin_target.py` amendment (C3 integration)

C3 (from [task 18](18_bin_target.md)) currently emits `target.fasta`
as the C3-selected primary contig for `plant_pt`. This task adds a
small hook on the `plant_pt` branch, placed after C3 has written
`target.fasta` and `secondaries.tsv` and before it writes
`bin_metadata.json`.

Pseudocode of the added section:

```
IF meta.assembly_target IS "plant_pt":
    import canonicalise_plastid from plastid_canonicalise   # local import
    iso_dir <- <out_root>/plastid_isoforms
    result  <- canonicalise_plastid(gfa_path, iso_dir)

    bin_metadata["plastid_canonicalisation"] <- fields of result
    bin_metadata["plastid_canonicalisation"]["substitution_applied"]
        <- (result.branch IS "canonical")

    IF result.branch IS "canonical":
        copy bytes of <iso_dir>/path1.fasta OVER <out_root>/target.fasta
    ELSE:
        # Safety net — C4 writes nothing off the canonical branch, so
        # an existing directory here means something went wrong.
        remove iso_dir if it exists (must be empty of isoform FASTAs)
ELSE:
    bin_metadata["plastid_canonicalisation"] <- {
        branch: "not_applicable",
        reason: "assembly_target=<meta.assembly_target>",
    }
```

Notes:

- **Import at call site, not module top.** Keeps C3's non-plant_pt
  code paths free of the C4 dependency (defensive — the container
  will always have C4 available, but the import is scope-local so a
  hypothetical C4-missing environment doesn't crash animal_mt runs).
- **Substitution is a byte copy, not a re-open.** `path1.fasta` is
  already valid single-record FASTA; copying bytes preserves the
  exact header written by C4. Downstream stages don't care about the
  record ID — they run on sequence.
- **`bin_metadata.json` enrichment is additive.** Task 18's existing
  metadata keys are untouched; only the `plastid_canonicalisation`
  sub-object is added. The key name is fixed by
  [spec §4.2](../../spec/plastid-canonicalisation.md).

## 3. `modules/local/bin_target.nf` — uncomment isoforms emit

Task 18 left the `emit: isoforms` block commented at
[modules/local/bin_target.nf:20-21](../../modules/local/bin_target.nf#L20-L21).
Uncomment it, so the output block reads:

```groovy
output:
tuple val(meta), path("target.fasta"), path("secondaries.tsv"), emit: binned
path("bin_metadata.json"), emit: metadata
path("plastid_isoforms/"), optional: true, emit: isoforms
```

The `optional: true` handles the animal_mt / plant_mt / non-canonical
cases where C3 does not create the directory.

No `script:` change — C3 invokes C4 in-process; the module just needs
to expose the resulting directory as an output.

## 4. Tighten `plant_pt` `bin_bounds.json`

[Task 18 §6](18_bin_target.md) set the plant_pt `target_min_bp` to
`70000` to accept LSC-only pre-canonicalisation. This task tightens
it now that `target.fasta` is guaranteed to be the ~150 kb
canonicalised path1 on the SRR-derived plastid fixture:

```json
{
    "target_min_bp": 140000,
    "target_max_bp": 170000,
    "target_max_contigs": 1
}
```

`target_max_contigs` tightens from 3 → 1 because path1 is a single
record. If the integration fixture happens to Flye-resolve as a
single circle (1-edge branch, no substitution), the resolved circle
is also single-record ~150 kb — the bound holds.

If the fixture resolves as N≠3 (non_canonical, no substitution), the
task will fail the assertion — that is a real biological signal
(unresolved plastid graph) that should be investigated rather than
silently loosened. Document the observation in Outcomes and file a
follow-up in `todo.md`.

## 5. Unit tests — `scripts/tests/test_plastid_canonicalise.py`

Independent of Nextflow. Fixtures constructed inline as strings — no
checked-in binary GFAs.

**Cases:** implement the full test-case matrix in
[spec §9](../../spec/plastid-canonicalisation.md) (3-edge canonical;
upper-case depth tag; missing depth tags; LSC–IR collision; 1-edge
resolved circle; 2-edge; 4-edge; the round-trip length invariant
`len(path1) == LSC + 2·IR + SSC`; and the `path2` reverse-complement
offset check). Assert the exact `branch` and `non_canonical_reason`
values the spec prescribes, and that no FASTA files are written off
the canonical branch.

**Fixture builder:** a helper that takes a list of
`(edge_name, length, depth)` tuples and renders a minimal GFA string.
Sequences are generated deterministically per edge — distinct
nucleotide content per edge so the sequence-identity assertions are
meaningful — so the same test always sees the same bytes. Keeps tests
fully self-contained.

**Coverage target:** 100 % branch coverage on the parsing,
classification and canonicalisation functions per
[CONSTITUTION.md rule 14](../../CONSTITUTION.md).

## 6. Integration-test wiring

Append a C4 block after the BIN_TARGET block in
[`tests/integration/assertions.sh`](../../tests/integration/assertions.sh),
following the existing block's style (read the branch from
`bin_metadata.json` with `jq`, set `FAILED=1` on failure, echo an
`OK:` / `FAIL:` / `WARN:` line per check).

Assertion logic, keyed on
`.plastid_canonicalisation.branch` for `INT-PLANT-01-pt`:

| Branch | Assertions |
|---|---|
| `canonical` | `plastid_isoforms/path1.fasta` and `path2.fasta` both exist and are non-empty; `target.fasta` is byte-identical to `path1.fasta` (`cmp -s`) |
| `resolved_circle` | `plastid_isoforms/` does **not** exist |
| `non_canonical` | `WARN` only, echoing the edge count — not a hard fail (see §11) |
| missing / unknown | `FAIL` — metadata absent or branch unrecognised |

Plus, for `INT-ANIMAL-01` and `INT-PLANT-01-mt`: `FAIL` if a
`plastid_isoforms/` directory exists at all.

Also:

- Tighten `tests/integration/expected/plant_pt/bin_bounds.json`
  per §4.
- Update the progressive-uncomment header in `assertions.sh` to note
  that C4 plastid canonicalisation is real.

**No new fixture files** beyond the tightened `bin_bounds.json`.

## 7. Fast CI

No change. C4 is invoked from within C3's execution path, but the stub
`bin_target.nf` script block only `touch`es outputs — it never
invokes `bin_target.py`. The `emit: isoforms` output is
`optional: true`, so its absence during `-stub-run` doesn't break
channel typing. `-stub-run` continues green.

## 8. Verification

1. `pytest scripts/tests/test_plastid_canonicalise.py` — every case
   in the spec matrix passes; coverage report shows 100 % branch on
   the classifier and canonicalisation functions.
2. `flake8 bin/plastid_canonicalise.py
   scripts/tests/test_plastid_canonicalise.py` — clean.
3. `nextflow run . -profile stub -stub-run` — 17/17 processes hit
   their stub, no regressions.
4. `nextflow run . -profile integration` locally:
   - `INT-PLANT-01-pt/bin_target/bin_metadata.json`
     `.plastid_canonicalisation.branch` is `canonical`.
   - `plastid_isoforms/path1.fasta` and `path2.fasta` both ~150 kb.
   - `target.fasta` byte-identical to `path1.fasta`.
   - `INT-ANIMAL-01` and `INT-PLANT-01-mt` have no
     `plastid_isoforms/` directory and
     `.plastid_canonicalisation.branch == "not_applicable"`.
5. `bash tests/integration/assertions.sh` — new C4 block reports OK
   for all three samples; the task 18 BIN_TARGET assertion also
   reports OK with the tightened plant_pt bounds.
6. CI green on `Tests` (fast) + nightly `Integration` on the PR.

## 9. Deliverables checklist

- [x] `bin/plastid_canonicalise.py` — original implementation written
      from
      [spec/plastid-canonicalisation.md](../../spec/plastid-canonicalisation.md)
      only; library + thin CLI.
- [x] `bin/bin_target.py` — plant_pt post-selection hook calling
      `canonicalise_plastid`; substitution + metadata enrichment.
      (Also added a `--gfa` CLI arg, needed to pass the assembly graph
      path through to C4 — not spelled out in the task pseudocode but
      required for the wiring to work; the `.nf` module already had
      `path(gfa)` as a process input from task 18.)
- [x] `scripts/tests/test_plastid_canonicalise.py` — full spec §9
      matrix (plus #2b for the `dp:i:` tag variant found during
      integration testing) with 100 % branch coverage on the
      classifier + canonicalisation functions (99% overall; the sole
      uncovered branch is the `if __name__ == "__main__":` CLI guard,
      outside the required scope).
- [x] [`modules/local/bin_target.nf`](../../modules/local/bin_target.nf) —
      `emit: isoforms` uncommented; `--gfa ${gfa}` added to the
      `script:` block; no other `script:` change.
- [x] [`tests/integration/expected/plant_pt/bin_bounds.json`](../../tests/integration/expected/plant_pt/bin_bounds.json) —
      `target_min_bp` tightened `70000 → 140000`;
      `target_max_contigs` tightened `3 → 1`.
- [x] [`tests/integration/assertions.sh`](../../tests/integration/assertions.sh) —
      C4 canonicalisation-branch block appended; progressive-uncomment
      header updated.
- [x] Fast CI equivalent (`nextflow run . -profile stub -stub-run`,
      17/17 processes) and `-profile integration` (47/47 processes,
      real tools/fixtures) both green locally. No PR opened in this
      session, so the actual GitHub Actions `Tests`/`Integration`
      workflows have not run against this branch yet.

## 10. Follow-up tasks

Append to [tasks/todo.md](../todo.md) on completion (planned as part of
the parent breakout plan; do not re-add if already present):

- When the `ORGANELLE_MAP` task is drafted, its rendering pass
  should walk both `path1` and `path2` when `plastid_isoforms/` is
  present (per [spec §3.6](../../spec/03-organelles.md)).

## 11. Notes / non-issues

- **No container change.** C4 sits in `bin/`; the C3 container built
  in task 18 already has Python 3.12 + `biopython`. C4 imports
  succeed because `bin/` is on `PATH` inside the container per
  Nextflow's staging convention.
- **Why `bin/` and not a package.** Keeping C4 as a flat file next to
  C3 matches the existing custom-logic layout
  ([bin/coverage_gate.py](../../bin/coverage_gate.py),
  [bin/parse_samplesheet.py](../../bin/parse_samplesheet.py)) — no
  package machinery, no `__init__.py`, C3 imports via a bare
  `from plastid_canonicalise import canonicalise_plastid`.
- **Why byte-copy, not in-place rename.** `path1.fasta` in
  `plastid_isoforms/` is the auditable artefact (recorded in
  `bin_metadata.json`); `target.fasta` at the root is the downstream
  contract. Keeping both files (with a byte-copy) makes provenance
  explicit — an auditor can `cmp` them and confirm the substitution.
- **Non-canonical is a WARN, not a FAIL.** A non-canonical plastid
  graph (edge count ≠ 1 or 3) is a real biological outcome — the
  operator needs to see it in the report, not have the pipeline
  hard-fail. The integration assertion treats it as an informative
  outcome to preserve
  [spec principle 7 (negative clarity)](../../CONSTITUTION.md).
- **Degenerate LSC == IR guard.** On a 3-edge input where the longest
  edge is also the deepest, the classifier returns `non_canonical`
  with the spec's collision reason rather than crashing on the
  remainder lookup. This never happens on real plastid GFAs (IR depth
  is ~2× the SC edges by construction), but the guard keeps the
  classifier total.

## 12. Outcomes

- Unit test cases pass: 25/25 (`scripts/tests/test_plastid_canonicalise.py`,
  run via `scripts/pytest.sh` inside `neoformit/daff-wf5-scripts:test`).
  Covers the full spec §9 matrix (#1–#10) plus #2b (see spec extension
  below), parser edge cases (duplicate names, `*` sequence, malformed
  float, non-first tag position), all four tie-break branches in
  `_select_lsc`/`_select_ir`, the zero-length-edge check, the
  construction-invariant guard (forced via a mocked `Seq`), and CLI
  smoke tests.
- Branch coverage on classifier + canonicalisation: 100 % on
  `parse_gfa`, `_select_lsc`, `_select_ir`, `_classify_three_edge`,
  `canonicalise_plastid` (99 % of the whole file; only the
  `if __name__ == "__main__":` guard is uncovered, which is outside
  the required scope).
- Observed canonicalisation branch on `INT-PLANT-01-pt`: `canonical`
  (edge_count=3), after the spec fix below. `path1_len == path2_len ==
  155277` bp, within the 145–160 kb expected range; `target.fasta`
  confirmed byte-identical to `plastid_isoforms/path1.fasta`
  (`cmp -s`).
- Assigned LSC / IR / SSC edge names on `INT-PLANT-01-pt`: `lsc_edge =
  edge_2`, `ir_edge = edge_1`, `ssc_edge = edge_3`.
- Any Flye 3-edge → non-canonical surprises: **Yes — a spec gap, not
  an algorithm bug.** The real Flye `assembly_graph.gfa` for
  `INT-PLANT-01-pt` tags depth as `dp:i:220` (integer GFA type), not
  `dp:f:56.7` (float) as spec/plastid-canonicalisation.md §3.3
  originally documented. The parser (correctly, per the original
  spec text) treated every edge's depth as absent → 0.0 → all three
  edges tied → `non_canonical` / `"depth_tie: all edges have equal
  depth"`. Per the task's clean-room instruction ("if something
  genuinely is not [in the spec], stop and file a follow-up to extend
  the spec rather than inferring the behaviour"), this was raised to
  the user rather than silently patched. With approval: extended
  spec/plastid-canonicalisation.md §3.3 (and §9, §10) to document that
  Flye may emit `dp:i:` as well as `dp:f:`, and that the parser must
  be type-character-agnostic on the GFA tag type; updated
  `DEPTH_TAG_RE` from `[Dd][Pp]:f:` to `[Dd][Pp]:[a-zA-Z]:`; added
  test #2b. Re-ran `-profile integration` from a clean `work/` after
  the fix and confirmed `canonical` (see above).
- Whether the tightened `140_000` lower bound needed adjustment: No —
  the real fixture's canonicalised path1 (155277 bp) sits comfortably
  inside `[140000, 170000]`.
- Two pre-existing task-18 issues, unrelated to C4, surfaced by this
  integration run and filed to [tasks/todo.md](../todo.md) rather than
  fixed here (out of this task's scope — task 18's `select_primary`
  / `check_circularity`, not C4): (1) `INT-ANIMAL-01`'s real assembly
  isn't detected as circular (`circular: false`); (2) both
  `INT-PLANT-01-pt` and `INT-PLANT-01-mt` have `n_target_selected ==
  0` — no contig clears C3's `target_candidate` thresholds on either
  real fixture. For `plant_pt` this is masked because C4's
  canonicalisation runs unconditionally on the raw
  `assembly_graph.gfa` and overwrites `target.fasta` on the
  `canonical` branch regardless of C3's own selection; `plant_mt` has
  no such rescue path, so its `assertions.sh` check fails
  (`target.fasta missing or empty`).
- Also filed to todo.md: `nextflow -resume` did not re-execute
  `BIN_TARGET` after an in-place edit to `bin/plastid_canonicalise.py`
  (the depth-regex fix) — had to `rm -rf work/ .nextflow*` and rerun
  from scratch to get a trustworthy result. Worth investigating
  whether this pipeline's cache mode actually hashes `bin/` script
  contents into the task hash.
