# Task 23 — Recalibrate `BIN_TARGET` (C3) against real-fixture evidence

**Phase:** P3 (from [spec §6](../../spec/06-phases.md)).
**Goal:** Bring [`bin/bin_target.py`](../../bin/bin_target.py) into line with
the revised binning criteria in
[spec §3.7](../../spec/03-organelles.md#37-target-binning-criteria-per-assembly-target),
which were written from the first `-profile integration` run to reach the
assertion stage (2026-07-31). Under the current criteria `plant_pt` and
`plant_mt` select **zero** contigs and `animal_mt` reports a circular
contig as linear.

**Prerequisite:** [task 18 — BIN_TARGET](18_bin_target.md) and
[task 20 — plastid canonicalisation](20_plastid_canonicalise.md)
are landed. This task changes C3's selection logic only; it does not
touch C4's algorithm.

**Not in scope:** the C4-substitution defect (see
[task 24](../24_plastid_substitution_guard.md)) and the COVERAGE_GATE
inflation defect (see [task 25](../25_coverage_gate_carryover.md)). Both were
found in the same investigation but are separate stages with separate
failure modes. Task 24 shares a file with this task — sequence 23 before
24.

## 1. Why the current criteria fail

Full evidence table in
[spec §3.7](../../spec/03-organelles.md#37-target-binning-criteria-per-assembly-target).
In summary, four defects, in descending order of severity:

| # | Defect | Consequence | Spec |
|---|---|---|---|
| 1 | No sibling-organelle discrimination | The two highest-coverage contigs in the `plant_mt` fixture are **chloroplast** (100 % of their length aligns to `plant_pt.mmi`). Ranking on coverage picks the plastid as the mitogenome. | [§3.7.1](../../spec/03-organelles.md#371-homology-is-measured-against-sibling-panels-too) |
| 2 | `aligned_frac` uses one alignment block | Understates contig coverage 2–10×; the sole reason `plant_pt` selects nothing. | [§3.7.2](../../spec/03-organelles.md#372-aligned-fraction-is-merged-across-all-alignment-blocks) |
| 3 | Coverage-spike gate assumes a nuclear background | Self-defeating downstream of RECRUIT, and **inverted** on `plant_mt`. | [§3.7.3](../../spec/03-organelles.md#373-coverage-is-a-ranking-signal-not-a-gate) |
| 4 | End-overlap circularity check | Structural false negative — Flye trims the overlap when it circularises. | [§3.7.4](../../spec/03-organelles.md#374-circularity-comes-from-flye-not-from-end-overlap) |

Defect 1 is the one that matters most: fixing 2–4 while leaving 1 in place
would emit 155 kb of chloroplast as a plant mitogenome — a confident,
plausible, wrong answer of exactly the kind
[rule 18](../../CONSTITUTION.md) exists to prevent.

## 2. Changes to `bin/bin_target.py`

### 2.1 Merged-interval homology metric

Replace `align_to_ref()`. For a contig and a loaded panel, collect **all**
`mappy` hits, merge their query intervals, and return merged aligned
fraction plus aggregate identity:

```
# pseudocode
hits         = list(aligner.map(seq))
merged       = merge_overlapping([(h.q_st, h.q_en) for h in hits])
aligned_frac = sum(e - s for s, e in merged) / len(seq)
identity     = 100 * sum(h.mlen for h in hits) / sum(h.blen for h in hits)
```

The interval merge is ~10 lines of stdlib and is the single most
test-worthy function in this task — give it direct unit tests
(adjacent, nested, disjoint, single, empty).

### 2.2 Sibling-panel discrimination

`BIN_TARGET` currently receives one `.mmi`. It needs the sibling
panel(s) per [§3.7.1](../../spec/03-organelles.md#371-homology-is-measured-against-sibling-panels-too).
The reference bundle already ships all three
([`refs/v<ver>/recruit/{animal_mt,plant_pt,plant_mt}.mmi`](../../spec/04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle)),
and [`modules/local/bin_target.nf`](../../modules/local/bin_target.nf) already
stages the bundle root and indexes into it — so this is a script-side
change with **no channel-topology change** and no new reference data.

Sibling map:

```
SIBLING_PANELS = {
    'plant_mt':  ['plant_pt'],
    'plant_pt':  ['plant_mt'],
    'animal_mt': [],
}
```

Record per-contig, in `secondaries.tsv` and `bin_metadata.json`, the
aligned fraction against **every** panel scored, not just the declared
one. An auditor looking at a rejected contig must be able to see *which*
organelle it was assigned to. New classification value:
`sibling_organelle` (distinct from `off_target`).

Guard the case where a sibling `.mmi` is absent from an older bundle:
warn to stderr, record `sibling_panels_scored: []` in metadata, and
continue — do not fail the sample
([principle 8](../../CONSTITUTION.md), cross-sample failure isolation).

### 2.3 Per-target thresholds

Replace the module-level `MIN_REF_IDENTITY` / `MIN_ALIGNED_FRAC` /
`COVERAGE_SPIKE_FACTOR` constants with the per-target table from
[§3.7.6](../../spec/03-organelles.md#376-revised-per-target-criteria).
Delete the coverage-spike admission gate; keep coverage as the
dominance sort key.

Per [rule 18](../../CONSTITUTION.md) and the standing to-do item on
hard-coded script config, these belong in `nextflow.config` as a
`params.bin_target_thresholds` map keyed by `assembly_target`, passed to
the script as CLI arguments — not as Python constants. This is the
natural point to make that move for C3, since the values are now
per-target and explicitly provisional pending
[§9 item 10](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking).

### 2.4 Circularity from Flye

Extend `parse_assembly_info()` to read the `circ.` column (field 4,
`Y`/`N`) — it currently stops at field 3. Then:

```
if flye_circ == 'Y':      circular, method = True,  'flye_circ'
elif check_circularity(): circular, method = True,  'end_overlap'
else:                     circular, method = False, 'none'
```

Record both `circular` and `circular_method` in `bin_metadata.json`.
Keep `check_circularity()` — it still covers contigs Flye failed to
circularise. Apply to all targets, not just `animal_mt`; circularity is
informative on `plant_pt` too, and the metadata field costs nothing.

### 2.5 Longest ORF demoted

Keep computing `orf_aa_len` and `selected_genetic_code` — both are
reported provenance ([rule 16](../../CONSTITUTION.md)) — but remove
`orf_ok` from the selection intersection per
[§3.7.5](../../spec/03-organelles.md#375-orf-integrity-as-implemented-is-vacuous).
Retain the column in `secondaries.tsv` as a diagnostic.

Removing a criterion without replacing it weakens C3's biological
evidence, which is why §3.7.5 names marker-gene presence as the
principled replacement. That replacement is carved out to
[task 26](../26_binning_marker_genes.md) — it needs a reference-bundle
addition and a [principle 6](../../CONSTITUTION.md) sign-off, neither of
which should block this task. **Do not implement marker-gene presence
here.**

## 3. Integration fixture corrections

`tests/integration/expected/plant_mt/bin_bounds.json` currently expects
200 000–700 000 bp. That bound was derived from the whole-assembly total,
which is now known to be **57 % plastid**. Hitting it requires binning
chloroplast as mitochondrion.

Correct it to the evidence: contig_3 + contig_5 ≈ 115 kb, so
`target_min_bp: 90000`, `target_max_bp: 400000`,
`target_max_contigs: 20`. Widen rather than pin tightly — the true plant
mt size for this fixture is not established (see §5.3).

Add to `tests/integration/assertions.sh`, in the BIN_TARGET block:

- every sample selects ≥ 1 contig (`n_target_selected >= 1`) — the
  assertion whose absence let `plant_pt` pass while selecting nothing;
- no emitted contig is classified `sibling_organelle`;
- `INT-PLANT-01-mt` does **not** emit contig_1 or contig_6 (the known
  plastid contigs) — a targeted regression assertion for defect 1;
- `INT-ANIMAL-01` reports `circular == true` **and**
  `circular_method == "flye_circ"`, so a future regression to
  end-overlap-only is caught rather than silently passing.

Update the progressive-uncomment header table in `assertions.sh` to
note the BIN_TARGET block was recalibrated by this task.

## 4. Unit tests — `scripts/tests/test_bin_target.py`

[Rule 14](../../CONSTITUTION.md): 100 % branch coverage of the selection
path. Existing cases 1–8 stay, with cases 3 and 4 re-baselined (the
coverage-spike gate they exercise is gone). New cases:

| # | Scenario | Expected |
|---|---|---|
| 9 | Interval merge: overlapping, adjacent, nested, disjoint, single, empty | Correct merged length in each case |
| 10 | Contig scoring higher against sibling panel than declared panel | `sibling_organelle`; excluded from `target.fasta` |
| 11 | Contig scoring higher against declared panel, both non-zero | `target_candidate` |
| 12 | Sibling `.mmi` missing from bundle | Warns, `sibling_panels_scored: []`, still selects, exit 0 |
| 13 | Flye `circ. = Y`, no end overlap | `circular true`, `circular_method "flye_circ"` |
| 14 | Flye `circ. = N`, synthetic 300 bp end overlap | `circular true`, `circular_method "end_overlap"` |
| 15 | Flye `circ. = N`, no overlap | `circular false`, `circular_method "none"` |
| 16 | `plant_mt` with 2 mt + 2 sibling-panel contigs | Exactly the 2 mt contigs emitted |
| 17 | All contigs below `min_aligned_frac` | Empty `target.fasta`, exit 0, all rows classified |
| 18 | Per-target thresholds applied from CLI args | Same contig admitted under `plant_mt` bounds, rejected under `animal_mt` bounds |

Case 16 is the regression test for the defect that motivates this task —
keep it synthetic (inline sequences) per
[rule 19](../../CONSTITUTION.md), no checked-in fixture data.

Run via `scripts/pytest.sh`; flake8 per the project venv.

## 5. Decisions required before / during implementation

These are the reason this task is written rather than the fix simply
applied. Each changes the shape of the work.

### 5.1 Marker-gene presence as criterion (c) — carved out

**Resolved: deferred to [task 26](../26_binning_marker_genes.md).** No
decision is required to land this task.

[§3.7.5](../../spec/03-organelles.md#375-orf-integrity-as-implemented-is-vacuous)
establishes the longest-ORF criterion is vacuous, and §2.5 removes it.
Replacing it with "contig carries target-appropriate organellar
protein-coding genes" is the biologically correct discriminator, but it
requires a broader binning gene set than the `plant_mt` barcode panel
(a reference-bundle addition, [rule 10](../../CONSTITUTION.md)), a C3
container rebuild ([rule 12](../../CONSTITUTION.md)), and an explicit
[principle 6](../../CONSTITUTION.md) sign-off. None of that should block
the four defect fixes in §2.

Task 26 opens with a go/no-go gate rather than an implementation: land
§2.1–2.5 first and confirm on real data how much of the problem the
sibling-panel test alone solves. The
[§3.7.6](../../spec/03-organelles.md#376-revised-per-target-criteria)
evidence suggests it may be sufficient — the `plant_mt` plastid contigs
score 1.000 against `plant_pt.mmi` versus 0.012–0.064 against
`plant_mt.mmi`. If selection is stable without marker genes, the bundle
change is not worth its maintenance cost
([rule 19](../../CONSTITUTION.md)) and task 26 closes unimplemented.

**In the meantime,** the NUMT/NUPT exposure this would have covered is
accepted and flagged — see §5.2.

### 5.2 NUMT/NUPT exposure

Recorded as accepted-and-flagged in
[§3.7.3](../../spec/03-organelles.md#373-coverage-is-a-ranking-signal-not-a-gate)
and tracked as
[§9 item 10a](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking).
Removing the coverage gate removes the incidental defence. For
`animal_mt` and `plant_pt` the dominance rule contains it (a NUMT will
not out-rank the real organelle on coverage). For `plant_mt`, which
emits all candidates, a NUMT clearing the homology floors **would** be
emitted.

Minimum mitigation for this task: where a `plant_mt` candidate's
coverage is below some fraction of the highest-ranked candidate's, emit
it but set a `low_coverage_candidate` flag in `secondaries.tsv` and
`bin_metadata.json` for the report to surface. Flagging, not filtering —
per [principle 7](../../CONSTITUTION.md), the operator sees the
uncertainty rather than the pipeline silently guessing. Confirm the
fraction (suggest 0.05) or confirm that flagging is sufficient.

### 5.3 Is the `plant_mt` fixture assembly actually complete?

The two genuine mt contigs total ~115 kb at 5–7× coverage, against a
§3.1 expectation of 200 kb – several Mb. The mitogenome is likely
under-assembled because the recruited read pool was dominated by
plastid — which is [task 25](../25_coverage_gate_carryover.md)'s subject.
Until that is resolved we cannot distinguish "C3 is missing contigs"
from "the assembly genuinely lacks them", which is why §3 widens the
bounds rather than pinning them. Do not tighten
`plant_mt/bin_bounds.json` until task 25 lands.

## 6. Exit criteria

- `-profile integration` selects ≥ 1 contig for all three fixtures;
  `INT-PLANT-01-mt` emits contig_3 + contig_5 and **not** contig_1 or
  contig_6.
- `INT-ANIMAL-01` reports `circular: true`,
  `circular_method: "flye_circ"`, and still selects contig_8.
- `INT-PLANT-01-pt` selects contig_2 **on C3's own evidence** —
  verifiable as `n_target_selected >= 1` in `bin_metadata.json`, not
  merely via C4's substitution (see [task 24](../24_plastid_substitution_guard.md)).
- `bin_metadata.json` records, per contig: aligned fraction against every
  panel scored, the winning panel, `circular_method`, and the thresholds
  applied.
- 100 % branch coverage on the selection path; `scripts/pytest.sh` green;
  flake8 clean at 79 cols.
- `-stub-run` still green — stub block unchanged.
- Thresholds live in `nextflow.config`, not in `bin/bin_target.py`.
- `spec/03-organelles.md` §3.7.6 updated if any threshold moved during
  implementation, with the observed value that motivated the move.

## 7. Follow-up

- Remove the three superseded `tasks/todo.md` entries (the
  `INT-ANIMAL-01` circularity item and the two `n_target_selected == 0`
  items) — this task and tasks 24–25 replace them.
- [§9 item 10 / 10a](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
  benchmarking remains open after this task; these are prototype
  defaults on three fixtures.
- [Task 26](../26_binning_marker_genes.md) — marker-gene presence as
  criterion (c). Its §2 go/no-go gate consumes this task's behaviour on
  the P3 datasets, so it cannot be scoped until this lands.

## Outcomes

Completed 2026-08-02. All of §2.1–2.5, §3 and §4 landed as specified.
No threshold in [§3.7.6](../../spec/03-organelles.md#376-revised-per-target-criteria)
required adjustment during implementation.

### Decisions taken

- **§5.2 NUMT/NUPT flag** — confirmed with the user at `0.05`.
  `low_coverage_fraction` is a per-target entry in
  `params.bin_target_thresholds`, and `low_coverage_candidate` is
  recorded in both `secondaries.tsv` and `bin_metadata.json`. No
  candidate in the current fixtures trips it (the two `plant_mt`
  candidates sit at 5× and 7×, a 0.71 ratio).
- **§5.3** — `plant_mt/bin_bounds.json` widened to 90 000–400 000 bp
  and deliberately **not** tightened, pending
  [task 25](../25_coverage_gate_carryover.md).

### Deviations from the brief

1. **`--organelle-ref` replaced by `--ref-dir`.** §2.2 called for
   sibling panels without a channel-topology change. C3 now takes the
   staged bundle directory and resolves `${panel}.mmi` itself, rather
   than taking one pre-resolved `.mmi` path. `modules/local/bin_target.nf`
   still stages the same `path organelle_refs` input — no channel change,
   as specified.
2. **`load_panel()` uses `seq_names`, not truthiness.** `mappy.Aligner`
   is truthy even for a missing, empty or corrupt index (verified in the
   C3 image), so the original `if not aligner` guard was unreachable
   dead code. An empty `aligner.seq_names` is the reliable load-failure
   signal; both the missing-file and corrupt-file paths are now covered
   by tests.
3. **`mappy` added to `scripts/requirements-test.txt`, and a compiler to
   the `TEST=1` Docker layer.** `scripts/tests/test_bin_target.py` could
   not previously run at all — the shared `neoformit/daff-wf5-scripts`
   image has no `mappy` (it ships in the separate BIN_TARGET image
   alongside minimap2), so importing `bin/bin_target.py` raised
   `ModuleNotFoundError`. Pinned to `mappy==2.31`, matching the C3
   container's minimap2 2.31. `mappy` is source-only on PyPI, so the
   test layer installs `gcc`/`zlib1g-dev`, builds, then purges them.
   **The runtime image (`TEST=0`) is untouched** and
   `requirements.txt` — the CI image-rebuild trigger — is unchanged.
4. **Case 2 re-baselined as well as cases 3 and 4.** §4 named 3 and 4,
   but case 2 ("only one passes ORF") tested a criterion §2.5 deletes.
   It now distinguishes a candidate from an `off_target` on homology.
5. **§7's three superseded `tasks/todo.md` entries were already
   removed** in commit `db2ae37` (the task-creation commit). No action
   needed.

### Verification

`-profile integration` (2026-08-02, `-resume`; BIN_TARGET and all
downstream stages re-executed). Every §6 exit criterion met:

| Sample | Selected | bp | Notes |
|---|---|---|---|
| `INT-ANIMAL-01` | contig_8 | 16 952 | `circular: true`, `circular_method: "flye_circ"` |
| `INT-PLANT-01-pt` | contig_2 | 155 277 | `n_target_selected: 1` — C3's own evidence, before C4 substitution |
| `INT-PLANT-01-mt` | contig_3 + contig_5 | 114 847 | contig_1 / contig_6 classified `sibling_organelle` |

The `plant_mt` per-panel fractions reproduce the
[§3.7](../../spec/03-organelles.md#37-target-binning-criteria-per-assembly-target)
evidence table exactly:

| Contig | cov | vs `plant_mt` | vs `plant_pt` | Classification |
|---|---|---|---|---|
| contig_3 | 7× | 0.3529 | 0.1212 | `target_candidate` |
| contig_5 | 5× | 0.1654 | 0.0872 | `target_candidate` |
| contig_1 | 168× | 0.0639 | **1.0000** | `sibling_organelle` |
| contig_6 | 137× | 0.0121 | **1.0000** | `sibling_organelle` |

- `bash tests/integration/assertions.sh` — all assertions passed,
  including the four added by §3.
- `-profile stub -stub-run` — green, stub block unchanged.
- `scripts/pytest.sh scripts/tests/test_bin_target.py` — 46 passed,
  98 % branch coverage of `bin/bin_target.py`. The three uncovered
  statements are the C4 `shutil.copyfile`/`rmtree` substitution
  branches ([task 24](../24_plastid_substitution_guard.md)'s scope) and the
  `__main__` guard; the selection path is 100 %.
- `flake8 bin/ scripts/` clean at 79 columns.

### Note for task 24 / follow-up

`scripts/tests/test_coverage_gate.py` fails (4 tests) in the shared
scripts test image because `seqkit` is not installed there. Pre-existing
and unrelated to this task, but it means `scripts/pytest.sh` is not
currently green as a whole suite — logged in `tasks/todo.md`.

`INT-PLANT-01-pt`'s `circular_method` is `end_overlap`, not
`flye_circ`: Flye marks contig_2 `N`, and the fallback fires on a
genuine residual overlap. That is the §3.7.4 fallback behaving as
designed and is recorded, not asserted.
