# Task 20 — P3 stage 10 helper: `plastid_canonicalise` implementation (C4)

**Phase:** P3 (from [spec §6](../spec/06-phases.md)).
**Goal:** Implement the plastid quadripartite canonicalisation
algorithm — described in
[spec/plastid-canonicalisation.md](../spec/plastid-canonicalisation.md)
(produced by [task 19](19_plastid_algorithm_spec.md)) — as
`bin/plastid_canonicalise.py` (custom-logic component C4 in
[spec §2.2](../spec/02-stages.md#22-custom-logic-components)); wire
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

## Clean-room firewall — read only the algorithm spec

[Task 19](19_plastid_algorithm_spec.md) produced
[spec/plastid-canonicalisation.md](../spec/plastid-canonicalisation.md)
by reading `reference-material/ptgaul/combine_gfa.py` and translating
the algorithm into a self-sufficient prose + pseudocode +
test-case-matrix document. **This task's implementer works solely
from that spec.**

Rules for the implementer of this task:

- **Do not open** `reference-material/ptgaul/combine_gfa.py`. If a
  detail is missing from the algorithm spec, file a follow-up to
  extend the spec (via task 19) rather than consulting the ptGAUL
  source.
- The finalised spec at
  [spec/plastid-canonicalisation.md](../spec/plastid-canonicalisation.md)
  is the only permitted algorithm reference. If task 19's Outcomes
  section references specific ptGAUL constructs, treat those as
  provenance metadata, not as an implementation reference.
- Do **not** import the ptGAUL script, vendor it under `bin/`, or
  paste any fragment (variable names, comment structure, function
  shape) from it into our source.
- Follow the CLI arg-surface requirements from
  [spec/plastid-canonicalisation.md](../spec/plastid-canonicalisation.md)
  — designed independently of ptGAUL's `-e / -d / -o` triple.

**Rationale:** two-task clean-room design is the industry-standard
defence against inadvertent code copying from unlicensed sources.
Task 19 does the reading; task 20 does the writing; the algorithm
spec is the firewall between them. Neither task alone touches both
sides of the wall.

---

**Prerequisites (both hard-blockers):**

1. [task 19 — plastid algorithm spec](19_plastid_algorithm_spec.md)
   is complete and its output
   ([spec/plastid-canonicalisation.md](../spec/plastid-canonicalisation.md))
   is committed. Without the spec, this task cannot start — the
   implementer has no permitted algorithm reference.
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
- `bin_metadata.json.canonicalisation` records edge count, branch
  (`canonical` / `resolved_circle` / `non_canonical`), the assigned
  LSC / IR / SSC edge IDs on `canonical`, and whether the substitution
  was applied.
- `animal_mt` and `plant_mt` samples produce **no**
  `plastid_isoforms/` directory (guard:
  `meta.assembly_target == 'plant_pt'`).
- `tests/integration/expected/plant_pt/bin_bounds.json`
  `target_min_bp` tightened from `70000` → `140000` (see §4).
- `bin/plastid_canonicalise.py` has unit tests covering all four GFA
  branches (3-edge canonical, 1-edge circle, 2-edge, 4+-edge) using
  synthetic GFAs, plus a round-trip length invariant
  (`len(path1) == LSC + 2·IR + SSC`), with 100 % branch coverage on
  the classifier and canonicalisation functions per
  [CONSTITUTION.md rule 14](../CONSTITUTION.md).
- Integration assertion added to
  [`tests/integration/assertions.sh`](../tests/integration/assertions.sh)
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
  [spec §3.6 step 5](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation-ptgaul-derived).
- Rendering both isoforms in `ORGANELLE_MAP` — deferred to the
  ORGANELLE_MAP task drafting, tracked in [tasks/todo.md](todo.md).
- Handling of GFA `L` (link) lines — canonicalisation reads only `S`
  lines (edge sequence + depth). Link orientation is BandageNG's
  concern.

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- No container change — inherit task 18's SHA pin.
- No host tools; Python + `biopython` come from the C3/C4 shared
  container.
- `bin/plastid_canonicalise.py` is auto-staged onto every process
  `PATH` by Nextflow's `bin/` convention — importable by C3 because
  both scripts sit under `bin/`.

---

## 1. `bin/plastid_canonicalise.py` (new — custom logic C4)

Python + `biopython` (already in the C3 container). Implementation
follows [spec/plastid-canonicalisation.md](../spec/plastid-canonicalisation.md);
if the pseudocode below and the spec disagree, the spec wins.

Pseudocode of the intended shape (illustrative — the algorithm
specification is authoritative):

```python
#!/usr/bin/env python3
"""C4 — plastid quadripartite canonicalisation.

Implements the algorithm described in
spec/plastid-canonicalisation.md. Parses a Flye assembly_graph.gfa,
classifies edge count, and (on the canonical 3-edge case) emits the
two SSC-orientation isoforms path1 and path2.

Callable as a library (canonicalise_gfa) from bin/bin_target.py, or
standalone from the CLI for debugging. Written from the algorithm
spec only; upstream ptGAUL source is not consulted or vendored.
"""

from Bio.Seq import Seq
from pathlib import Path
from typing import NamedTuple
import argparse
import json
import re
import sys


class Edge(NamedTuple):
    edge_id: str
    sequence: Seq
    depth: float
    length: int


class Result(NamedTuple):
    branch: str          # 'canonical' | 'resolved_circle' | 'non_canonical'
    edge_count: int
    lsc_id: str | None
    ir_id: str | None
    ssc_id: str | None
    path1_len: int | None
    path2_len: int | None


# --- GFA parsing -----------------------------------------------------

def parse_edges(gfa_path: Path) -> list[Edge]:
    """Read S lines; extract sequence and depth tag (dp:f:N | DP:f:N)."""
    edges: list[Edge] = []
    for line in gfa_path.read_text().splitlines():
        if not line.startswith("S"):
            continue
        parts = line.split("\t")
        edge_id, seq_str = parts[1], parts[2]
        depth = _extract_depth(parts[3:])
        edges.append(Edge(edge_id, Seq(seq_str), depth, len(seq_str)))
    return edges


def _extract_depth(tags: list[str]) -> float:
    """Return float value of the first dp:f:N or DP:f:N tag; 0 if absent."""
    for tag in tags:
        m = re.match(r"^(?:dp|DP):f:([-\d.eE+]+)$", tag)
        if m:
            return float(m.group(1))
    return 0.0


# --- Classification + canonicalisation -------------------------------

def canonicalise_gfa(gfa_path: Path, outdir: Path) -> Result:
    """Public library entry point invoked from bin/bin_target.py.

    Writes path1.fasta and path2.fasta into outdir on the canonical
    branch; otherwise writes nothing. Returns a Result summary either
    way.
    """
    edges = parse_edges(gfa_path)
    n = len(edges)
    if n == 1:
        return Result("resolved_circle", 1, None, None, None, None, None)
    if n != 3:
        return Result("non_canonical", n, None, None, None, None, None)

    lsc = max(edges, key=lambda e: e.length)
    ir = max(edges, key=lambda e: e.depth)
    if ir.edge_id == lsc.edge_id:
        # Degenerate — longest edge is also deepest.
        return Result("non_canonical", n, None, None, None, None, None)
    ssc = next(e for e in edges if e.edge_id not in {lsc.edge_id, ir.edge_id})

    ir_rc = ir.sequence.reverse_complement()
    ssc_rc = ssc.sequence.reverse_complement()

    path1 = lsc.sequence + ir.sequence + ssc.sequence + ir_rc
    path2 = lsc.sequence + ir.sequence + ssc_rc + ir_rc

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "path1.fasta").write_text(f">path1\n{path1}\n")
    (outdir / "path2.fasta").write_text(f">path2\n{path2}\n")

    return Result(
        "canonical", n,
        lsc.edge_id, ir.edge_id, ssc.edge_id,
        len(path1), len(path2),
    )


# --- CLI (independent of upstream ptGAUL arg shape) ------------------

def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("gfa", type=Path)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--json-out", type=Path)
    args = p.parse_args()

    result = canonicalise_gfa(args.gfa, args.outdir)
    if args.json_out:
        args.json_out.write_text(json.dumps(result._asdict(), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

Design notes:

- **`biopython` only.** No new deps — the C3 container from task 18
  already ships `biopython` for ORF search.
- **Pure library + thin CLI.** `canonicalise_gfa()` is the tested
  interface; the CLI is a wrapper for standalone debugging. C3
  imports the function directly.
- **Depth tag regex.** Flye writes `dp:f:<value>` on `S` lines; the
  fallback `DP:f:` covers assemblers that use the upper-case
  convention. Missing tag → depth 0 (won't be selected as IR).
- **Degenerate-input guard.** If the longest and deepest edges are
  the same (impossible on a real plastid, but common on synthetic
  1-edge or misclassified graphs), fall to `non_canonical` rather
  than crash on a `next(...)` StopIteration.

## 2. `bin/bin_target.py` amendment (C3 integration)

C3 (from [task 18](18_bin_target.md)) currently emits
`target.fasta` as the C3-selected primary contig for `plant_pt`. This
task adds a small hook on the `plant_pt` branch:

Pseudocode of the added section:

```python
# In bin_target.py, after C3 has written target.fasta and
# secondaries.tsv, and before writing bin_metadata.json:

if meta_assembly_target == "plant_pt":
    from plastid_canonicalise import canonicalise_gfa
    iso_dir = out_root / "plastid_isoforms"
    canon = canonicalise_gfa(gfa_path, iso_dir)
    bin_metadata["canonicalisation"] = canon._asdict()
    bin_metadata["canonicalisation"]["substitution_applied"] = (
        canon.branch == "canonical"
    )
    if canon.branch == "canonical":
        # Overwrite target.fasta with path1 content.
        (out_root / "target.fasta").write_bytes(
            (iso_dir / "path1.fasta").read_bytes()
        )
    else:
        # No isoforms directory should exist on non-canonical branches.
        if iso_dir.exists():
            # Only path1/path2 could be there; canonicalise_gfa never
            # writes on non-canonical, so this branch is a safety net.
            for f in iso_dir.iterdir():
                f.unlink()
            iso_dir.rmdir()
else:
    bin_metadata["canonicalisation"] = {
        "branch": "not_applicable",
        "reason": f"assembly_target={meta_assembly_target}",
    }
```

Notes:

- **Import at call site, not module top.** Keeps C3's non-plant_pt
  code paths free of the C4 dependency (defensive — the container
  will always have C4 available, but the import is scope-local so a
  hypothetical C4-missing environment doesn't crash animal_mt runs).
- **Substitution is a byte copy, not a re-open.** `path1.fasta` is
  already valid single-record FASTA; copying bytes preserves the
  exact `>path1\n<seq>\n` header. Downstream stages don't care about
  the record ID — they run on sequence.
- **`bin_metadata.json` enrichment is additive.** Task 18's
  existing metadata keys are untouched; only the
  `canonicalisation` sub-object is added.

## 3. `modules/local/bin_target.nf` — uncomment isoforms emit

Task 18 left the `emit: isoforms` block commented at
[modules/local/bin_target.nf:20-21](../modules/local/bin_target.nf#L20-L21).
Uncomment it:

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

Independent of Nextflow. Fixtures constructed inline as strings —
no checked-in binary GFAs.

Cases:

| # | Input | Expected branch | Extra assertions |
|---|---|---|---|
| 1 | 3-edge synthetic: LSC 90 kb (depth 1), IR 25 kb (depth 2), SSC 15 kb (depth 1) | `canonical` | `path1.fasta` + `path2.fasta` written; both length 155 kb; `result.lsc_id/ir_id/ssc_id` correct |
| 2 | Same as 1 but with `DP:f:` (upper-case) tag | `canonical` | Depth still parsed, IR correctly identified |
| 3 | Same as 1 but IR missing depth tag | `non_canonical` (depth 0 makes LSC also "deepest", degenerate guard fires) | No files written |
| 4 | 1-edge synthetic (`resolved_circle`) | `resolved_circle` | No files written |
| 5 | 2-edge synthetic | `non_canonical` | No files written |
| 6 | 4-edge synthetic | `non_canonical` | No files written |
| 7 | 3-edge, LSC and IR happen to be same edge (synthetic degenerate) | `non_canonical` | No files written; safety-net guard exercised |
| 8 | Round-trip: `len(path1) == LSC.length + 2*IR.length + SSC.length` | canonical | Case 1's `result.path1_len == 155000` |
| 9 | `path2` uses reverse-complement of SSC (not SSC itself) | canonical | Read `path2.fasta`, extract SSC region (chars LSC_len+IR_len : LSC_len+IR_len+SSC_len), verify `== rc(SSC.sequence)` |

Fixture builder: a helper `make_gfa(edges)` where `edges` is a list of
`(edge_id, length, depth, seed)` and sequences are generated
deterministically from `seed` so the same test always sees the same
bytes. Keeps tests fully self-contained.

**Coverage target:** 100 % branch coverage on `canonicalise_gfa`,
`parse_edges`, and `_extract_depth` per
[CONSTITUTION.md rule 14](../CONSTITUTION.md).

## 6. Integration-test wiring

Append after the BIN_TARGET block in
[`tests/integration/assertions.sh`](../tests/integration/assertions.sh):

```bash
# C4 plastid canonicalisation is real (task 20) — plant_pt only:
sample="INT-PLANT-01-pt"
iso_dir="$OUTDIR/$sample/bin_target/plastid_isoforms"
path1="$iso_dir/path1.fasta"
path2="$iso_dir/path2.fasta"
tgt="$OUTDIR/$sample/bin_target/target.fasta"
meta="$OUTDIR/$sample/bin_target/bin_metadata.json"

branch=$(jq -r .canonicalisation.branch "$meta" 2>/dev/null || echo "missing")
case "$branch" in
    canonical)
        for f in "$path1" "$path2"; do
            if [[ ! -s "$f" ]]; then
                echo "FAIL: $sample canonical branch but $f missing"
                FAILED=1
            fi
        done
        # target.fasta must be byte-identical to path1
        if ! cmp -s "$tgt" "$path1"; then
            echo "FAIL: $sample target.fasta not byte-identical to path1.fasta"
            FAILED=1
        else
            echo "OK:   $sample plastid canonicalisation applied (target == path1)"
        fi
        ;;
    resolved_circle)
        if [[ -d "$iso_dir" ]]; then
            echo "FAIL: $sample resolved_circle branch but plastid_isoforms/ exists"
            FAILED=1
        else
            echo "OK:   $sample plastid resolved as single circle (no isoforms)"
        fi
        ;;
    non_canonical)
        echo "WARN: $sample plastid non-canonical (edge count $(jq .canonicalisation.edge_count "$meta")) — investigate"
        # Not a hard fail — this is real biological signal, not a code bug.
        ;;
    *)
        echo "FAIL: $sample canonicalisation branch unknown or metadata missing"
        FAILED=1
        ;;
esac

# animal_mt / plant_mt must NOT have a plastid_isoforms directory
for sample in INT-ANIMAL-01 INT-PLANT-01-mt; do
    if [[ -d "$OUTDIR/$sample/bin_target/plastid_isoforms" ]]; then
        echo "FAIL: $sample has plastid_isoforms/ but is not plant_pt"
        FAILED=1
    fi
done
```

Also:

- Tighten `tests/integration/expected/plant_pt/bin_bounds.json`
  per §4.
- Update the progressive-uncomment header in `assertions.sh` to note
  that C4 plastid canonicalisation is real.

**No new fixture files** beyond the tightened `bin_bounds.json`.

## 7. Fast CI

No change. C4 is invoked from within C3's stub-time execution path,
but the stub `bin_target.nf` script block only `touch`es outputs —
it never invokes `bin_target.py`. The `emit: isoforms` output is
`optional: true`, so its absence during `-stub-run` doesn't break
channel typing. `-stub-run` continues green.

## 8. Verification

1. `pytest scripts/tests/test_plastid_canonicalise.py` — all 9 cases
   pass; coverage report shows 100 % branch on the classifier and
   canonicalisation functions.
2. `flake8 bin/plastid_canonicalise.py
       scripts/tests/test_plastid_canonicalise.py
       --max-line-length=100` — clean.
3. `nextflow run . -profile stub -stub-run` — 17/17 processes hit
   their stub, no regressions.
4. `nextflow run . -profile integration` locally:
   - `INT-PLANT-01-pt/bin_target/bin_metadata.json.canonicalisation.branch`
     is `canonical`.
   - `plastid_isoforms/path1.fasta` and `path2.fasta` both ~150 kb.
   - `target.fasta` byte-identical to `path1.fasta`.
   - `INT-ANIMAL-01` and `INT-PLANT-01-mt` have no
     `plastid_isoforms/` directory and
     `canonicalisation.branch == "not_applicable"`.
5. `bash tests/integration/assertions.sh` — new C4 block reports OK
   for all three samples; the task 18 BIN_TARGET assertion also
   reports OK with the tightened plant_pt bounds.
6. CI green on `Tests` (fast) + nightly `Integration` on the PR.

## 9. Deliverables checklist

- [ ] `bin/plastid_canonicalise.py` — implementation written from
      [spec/plastid-canonicalisation.md](../spec/plastid-canonicalisation.md)
      only (no ptGAUL source or task 19 draft notes referenced during
      authoring); library + thin CLI.
- [ ] `bin/bin_target.py` — plant_pt post-selection hook calling
      `canonicalise_gfa`; substitution + metadata enrichment.
- [ ] `scripts/tests/test_plastid_canonicalise.py` — 9 cases with
      100 % branch coverage on the classifier + canonicalisation
      functions.
- [ ] [`modules/local/bin_target.nf`](../modules/local/bin_target.nf) —
      `emit: isoforms` uncommented; no `script:` change.
- [ ] [`tests/integration/expected/plant_pt/bin_bounds.json`](../tests/integration/expected/plant_pt/bin_bounds.json) —
      `target_min_bp` tightened `70000 → 140000`;
      `target_max_contigs` tightened `3 → 1`.
- [ ] [`tests/integration/assertions.sh`](../tests/integration/assertions.sh) —
      C4 canonicalisation-branch block appended; progressive-uncomment
      header updated.
- [ ] Fast CI (`Tests`) + `Integration` workflows both green on the PR.

## 10. Follow-up tasks

Append to [tasks/todo.md](todo.md) on completion (planned as part of
the parent breakout plan; do not re-add if already present):

- When the `ORGANELLE_MAP` task is drafted, its rendering pass
  should walk both `path1` and `path2` when `plastid_isoforms/` is
  present (per
  [spec §3.6 step 5](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation-ptgaul-derived)).

## 11. Notes / non-issues

- **No container change.** C4 sits in `bin/`; the C3 container built
  in task 18 already has Python 3.12 + `biopython`. C4 imports
  succeed because `bin/` is on `PATH` inside the container per
  Nextflow's staging convention.
- **Why `bin/` and not a package.** Keeping C4 as a flat file next to
  C3 matches the existing custom-logic layout
  ([bin/coverage_gate.py](../bin/coverage_gate.py),
  [bin/parse_samplesheet.py](../bin/parse_samplesheet.py)) — no
  package machinery, no `__init__.py`, C3 imports via bare
  `from plastid_canonicalise import canonicalise_gfa`.
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
  [spec principle 7 (negative clarity)](../CONSTITUTION.md).
- **Degenerate LSC == IR guard.** On a synthetic 3-edge input where
  the longest edge is also the deepest, `canonicalise_gfa` returns
  `non_canonical` rather than crash on the `next(...)` remainder
  lookup. This never happens on real plastid GFAs (IR depth is ~2×
  the SC edges by construction), but the guard keeps the classifier
  total.

## 12. Outcomes

- Unit test cases pass: —
- Branch coverage on classifier + canonicalisation: —
- Observed canonicalisation branch on `INT-PLANT-01-pt`: —
- Assigned LSC / IR / SSC edge IDs: —
- Any Flye 3-edge → non-canonical surprises: —
- Whether the tightened `140_000` lower bound needed adjustment: —
