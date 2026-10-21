# Task 64 — Protein panel: catch mislabelled RefSeq records at build time

**Phase:** P4b (plant annotation arms). The panel check applies to all
three origins. Rolled out of task 32_research_plastid_noncds.md §7.5.

**Status:** open. **No preconditions.** Independent of task
63_cds_locus_clustering.md and task 65_plastid_noncds_barrnap.md.
Should land before task 66_plastid_map_structure.md, whose visual
check would otherwise draw the phantom genes described below.

## Overview

The protein panel (`refs/<ver>/proteins/<origin>/<gene>.faa`, spec
§4.3) is what miniprot aligns against the assembly to find every
protein-coding gene. It is built from the CDS translations in RefSeq
organelle GenBank records. Each sequence's gene name is taken as-is
from the record's `/gene` qualifier, and nothing checks it.

Some RefSeq records are mislabelled. Two are confirmed in the current
`plant_pt` panel:

| Panel entry | Labelled as | Actually | Evidence |
|---|---|---|---|
| `NC_034951.1_ycf3` (*Solanum berthaultii*) | `ycf3` | `clpP` | 198 aa (ycf3 is ~168, clpP ~196). Starts `MPIGVPKV…`, the clpP N-terminus, not `MPRSRING…` like every other `ycf3`. Aligns exon-for-exon to the fixture's `clpP` locus. |
| `NC_058879.1_cemA` (*Camellia gigantocarpa*) | `cemA` | `petA` | 320 aa (cemA ~230, petA ~320). Starts `MQTRNTFSWIKEQ…`, matching `petA` (*S. lycopersicum* `MQTRNAFSWLKKQ…`). Calls a second "`cemA`" at the fixture's `petA` locus. |

A mislabelled record puts a **phantom gene** in the annotation: a
second, wrongly named feature on top of a real gene. The CDS merge
clusters by gene name, so it can't notice, and the organelle map
draws both. A mislabel on a barcode gene would be worse: a barcode
called from the wrong locus would go to Taxodactyl under the wrong
name. That is exactly the "subtly wrong answer that ships as result"
that rule 18 forbids.

The fix belongs in the panel build, where the bad data comes in. Each
candidate representative is compared with the others. One that looks
more like a different gene than like the gene it claims to be is
excluded, and the exclusion is recorded. A new reference bundle,
`v2026.10`, ships the cleaned panel.

## 1. Constraints

- **Principle 10.** The bundle is immutable and versioned. The fix is
  a new bundle (`v2026.10`, refs/VERSIONING.md), never an edit to
  `v2026.09_1`. Every exclusion goes in `proteins/provenance.json`
  and from there into `manifest.json`.
- **Rule 18.** Exclusions are recorded, never silent. A gene left
  below the ≥3-candidate floor by exclusions drops out of the panel
  the way thin genes already do (spec §4.3), and that is recorded
  too.
- **Rule 12/19.** The pairwise comparison runs in a container already
  pinned for this project (the BLAST 2.17.0 image used by
  `BLAST_VALIDATE` and `scripts/refdata/build_validate.sh`). No new
  image, and no new Python dependency in `build_proteins.py`, which
  stays stdlib-only and host-side.
- **Deterministic.** The same RefSeq release and arguments give the
  same panel. Representative selection (task 48's breadth-first
  traversal) stays deterministic once flagged accessions are removed
  from the candidate pool.
- **Panel-competition check (spec §4.3).** Any panel change must not
  lose or degrade a barcode locus on any fixture.

## 2. Scope

### 2.1 The label check

Run once per origin, over the **selected** representatives of every
gene. That is ~800 sequences for `plant_pt`, ~400 for `plant_mt`, and
~14,300 for `animal_mt` (task 48: 1,100 per gene × 13).

1. Concatenate the origin's selected representatives, each tagged
   with its labelled gene.
2. Run an all-vs-all `blastp` inside the BLAST image (`-outfmt 6` with
   `bitscore`). Self-hits are ignored.
3. For each representative, find the best bitscore to any *other*
   representative of its **own** gene (`own`) and the best bitscore to
   a representative of any **other** gene (`other`, plus which gene).
4. **Flag** it when `other > own`, or when it has no own-gene hit at
   all and does have an other-gene hit.
5. Remove flagged accessions from that gene's candidate pool and
   **re-run selection**, so the budget is refilled from the next
   candidates in traversal order. Re-check the replacements. Repeat
   until a pass flags nothing. Cap the iterations (3 is plenty) and
   fail the build loudly if the cap is reached.

A representative with no hits at all to its own gene and none to any
other gene is not flagged. That is a divergent but plausible sequence,
not evidence of a mislabel. Count these in provenance.

**Genuine paralogs.** Some organelle genes are homologous to each
other: `psaA`/`psaB`, `psbC`/`psbD`, and `nad`/`ndh` subunits across
targets (cross-target pairs don't matter, since the check runs per
origin). A real `psaA` should still score higher against other `psaA`
than against `psaB`. If the check flags many members of a paralogous
pair rather than isolated outliers, stop and report. Do not loosen
the rule to make the count go away.

### 2.2 Where it runs

`build_proteins.py full` currently selects and writes the panel in one
step (spec §4.3, task 48). Split it so the check fits in between,
without moving BLAST into a stdlib script:

```
build_proteins.py full --emit-selection sel.json   # select, write nothing final
docker run blast ... all-vs-all → hits.tsv          # inside the pinned BLAST image
build_proteins.py label-check sel.json hits.tsv     # flag; if any, write exclusions
    → re-run select excluding flagged; loop until clean
build_proteins.py full --exclude exclusions.json    # final write + provenance
```

That's a sketch. Any shape works that keeps BLAST in the container,
the selection and flagging logic in tested stdlib Python, and the loop
deterministic. Document the exact build sequence in the module
docstring and in spec §4.3, since there is no orchestration script
for this build today (it is run by hand, per task 48).

### 2.3 Provenance

`proteins/provenance.json` gains, per origin, a `label_check` block:

```
"label_check": {
  "method": "all-vs-all blastp 2.17.0 within origin; flag if best other-gene bitscore > best own-gene bitscore",
  "iterations": 2,
  "excluded": [
    {"accession": "NC_034951.1", "labelled_gene": "ycf3",
     "best_other_gene": "clpP", "own_bitscore": 41.2, "other_bitscore": 389.0},
    ...
  ],
  "no_hits": 3
}
```

The bitscore values above are made up. `build_manifest.py` folds the
block into `manifest.json` the same way it folds the existing strata
block.

### 2.4 Bundle `v2026.10`

Follow the procedure task 48 used for `v2026.09_1` (its Outcomes
record each step):

- Rebuild `proteins/` only. Copy `recruit/`, `validate/` and
  `annotate/` from `v2026.09_1` unchanged, and say in the manifest
  that they were carried over.
- Write `scripts/refs-v2026.10.sha256`.
- Bump the bundle version in `conf/integration.config`,
  `.github/workflows/integration.yml`, and the usage examples in
  `deploy/azure/run-wf5-azure.sh` and `scripts/publish_refs_extracted.sh`.
- **Uploading to Azure** (`daffstandard/refdata-wf5/v2026.10/`, spec
  §4a) publishes the bundle. **Confirm with the user before
  uploading.** Then verify with `bash scripts/fetch_refs.sh v2026.10`
  on a clean scratch checkout.
- Do not retire `v2026.09_1` from Azure while any branch references
  it (refs/VERSIONING.md).

## 3. Testing

### 3.1 Unit tests (`scripts/tests/test_build_proteins.py`)

Keep the existing coverage level. Add tests with small hand-written
hit tables. No BLAST in the unit tests: the boundary is the hits file
(task 27).

- A representative whose best hit is another gene is flagged, with
  the other gene and both scores recorded.
- A representative whose own-gene and other-gene scores are equal is
  **not** flagged. The rule is strict `>`. Test it.
- No own-gene hits plus an other-gene hit is flagged. No hits at all
  is counted, not flagged.
- Self-hits are ignored.
- After exclusion, re-selection refills the budget from the next
  candidate. The replacement is the one the traversal would have
  picked next, and the result does not depend on input order.
- A gene pushed below the minimum by exclusions is dropped and
  recorded.
- The loop stops when a pass is clean, and fails at the iteration
  cap.
- The `label_check` provenance block has the expected shape and is
  folded into the manifest (`build_manifest.py` test, if one exists;
  add one if not).

### 3.2 Build-level acceptance

- **Both known mislabels are excluded:** `NC_034951.1` from `ycf3`
  and `NC_058879.1` from `cemA`.
- List every other exclusion, per origin, in Outcomes. For each, give
  the labelled and inferred gene and a one-line sanity check (length,
  N-terminus, or source record). If the list for any origin is longer
  than ~2 % of its representatives, stop and review before shipping.
  That many would point at a flaw in the rule, not at RefSeq.

### 3.3 Integration and panel competition

Run `-profile integration` on `v2026.09_1` and on `v2026.10`. Record:

- **Barcodes (contractual):** per fixture, every barcode locus still
  passes with the same coordinates, or a shift is explained by
  improved identity and coverage (spec §4.3 panel-competition rule).
- **`INT-PLANT-01-pt` annotation:** exactly one `ycf3` and one `cemA`
  feature (today two each). No other gene's copy count changes unless
  explained.
- **`INT-ANIMAL-01`:** annotation unchanged unless exclusions touched
  `animal_mt`. If they did, explain every difference.

Extend `tests/integration/expected/plant_pt/annotation_bounds.json`
(created by task 63. If task 63 hasn't landed, create it in the same
shape) with a `single_copy_genes` list: genes that must appear
**exactly once**, `ycf3` and `cemA`. Assert it in
`tests/integration/assertions.sh`. Both are single-copy LSC genes in
land-plant plastomes, so this is biology, not a fixture-tuned bound.

### 3.4 Stub run

Unaffected (the stub profile does not read the panel). Confirm it
passes.

## 4. Spec updates

- spec/04-reference-data.md §4.3: add a "Label check" bullet after
  "Per-gene curation" with the rule, where it runs, and the recorded
  provenance. Add the label check to the "Version pinning" list.
- spec/04-reference-data.md §4.4 and any bundle-version references:
  `v2026.10`.

## 5. Out of scope

- A runtime guard against two gene names on one locus (see task 63
  §2.4).
- Re-tuning panel breadth or strata (task 48, spec §9 item 11).
- Alias-table changes in `assets/organelle_gene_sets.json`, unless a
  flagged record turns out to be a legitimate alias rather than a
  mislabel. Record it if so.

## 6. Checklist

- [x] §2.1–2.2 label check implemented and loop deterministic
- [x] §2.3 provenance block written and folded into manifest
- [x] §3.1 unit tests; flake8 clean
- [x] §3.2 both known mislabels excluded; full exclusion list reviewed
- [x] §2.4 `v2026.10_1` built, checksummed; versions bumped
- [x] §2.4 upload confirmed with user, uploaded, fetch verified
- [x] §3.3 before/after recorded; no barcode lost or degraded;
      `single_copy_genes` asserted
- [x] §4 spec updated
- [ ] Outcomes section written

## Outcomes

**Status: complete.**
Work was done in git worktree `../wf5-task64` (branch
`task-64-label-check`) because another session was editing the main
tree (task 76).

### Deviations (approved)

- **Bundle name `v2026.10_1`, not `v2026.10`.** Task 68 had already
  published `v2026.10` (masked `animal_mt` recruit panel). Per
  refs/VERSIONING.md the new bundle is `v2026.10_1`, with `recruit/`,
  `validate/`, `annotate/` carried over from `v2026.10` (not
  `v2026.09_1`), so task 68's masking is kept. Only `proteins/` is
  rebuilt. Version references bumped in `conf/integration.config`,
  `.github/workflows/integration.yml`, `deploy/azure/run-wf5-azure.sh`,
  `scripts/publish_refs_extracted.sh`, `docs/azure/README.md`. Fixture
  version stays `v2026.10`. Checksum: `scripts/refs-v2026.10_1.sha256`.
  The manifest `getorganelle` release field is carried over as `0.0.1`.

### Implementation

- `build_proteins.py full` gains `--emit-selection DIR` and
  `--exclude EXCLUSIONS_JSON`; new `label_check` mode (exit 0 clean,
  2 flagged, 1 at the pass cap). `scripts/refdata/build_proteins_checked.sh`
  runs the loop (BLAST in the pinned 2.17.0 image, `-evalue 1e-3
  -max_target_seqs 500`). Sequence documented in the module docstring
  and spec §4.3.
- Exclusions are keyed on (origin, gene, accession), as the brief says.
- A gene dropped by exclusions is recorded in `dropped_genes`; none
  were dropped in this build. A gene already thin before exclusions
  still fails the build, as before.
- Tests: 11 new in `test_build_proteins.py`, 1 in new
  `test_build_manifest.py`; full suite 100 % branch coverage on `bin/`;
  flake8 clean.

### Build (RefSeq release 237, regenerated from the release gbff)

CDS counts match task 48 (1,309,988 plastid; 25,978 plant_mt; 205,087
animal_mt). Passes to clean: `animal_mt` 2, `plant_mt` 2, `plant_pt` 3.
No-hit representatives (counted, not flagged): animal_mt 100, plant_mt
5, plant_pt 4. Exclusions: plant_pt 9 / 800 (1.1 %), plant_mt 1 / 410,
animal_mt 1 / 14,300 — all under the 2 % review threshold.

| Origin | Accession | Labelled → inferred | Check |
|---|---|---|---|
| plant_pt | NC_034951.1 | ycf3 → clpP | 198 aa `MPIGVPKV…` (known mislabel) |
| plant_pt | NC_058879.1 | cemA → petA | 320 aa `MQTRNTFS…` (known mislabel) |
| plant_pt | NC_058880.1 | cemA → petA | same 320 aa petA-like record, Camellia semiserrata |
| plant_pt | NC_040106.1 | rpl32 → rpl23 | 93 aa `MDGIKYAV…`, identical to NC_029371.1 rpl23 |
| plant_pt | NC_058854.1 | rps19 → rpl2 | 159 aa vs median 92 (the 49 aa record is also present) |
| plant_pt | NC_020019.1, NC_022264.1 | petA → cemA | **false positives, see below** |
| plant_pt | NC_061388.1 | clpP → ycf3 | **false positive, see below** |
| plant_pt | NC_029371.1 | rpl23 → rpl32 | **false positive, see below** |
| plant_mt | NC_007886.1 | mttB → atp8 | not independently verified; phantom `orfB` at contig_6 removed from the INT-PLANT-01-mt annotation |
| animal_mt | NC_070094.1 | ATP8 → COX2 | 74 aa vs median 54; best hit COX2 (bitscore 70.9) |

**Victims.** Four exclusions are real genes flagged only because a
genuinely mislabelled representative scored higher against them than
their own-gene peers did (e.g. true `petA` scores 629 against other
petA but 649 against the mislabelled `cemA`, which is itself a petA).
The rule as specified flags every such record in the same pass, so
the mislabelled record and its victims go together. Cost: a good
representative is swapped for the next candidate in traversal order.
The provenance "excluded" lists therefore overstate mislabels by about
4. No paralog-pair mass flagging occurred (psaA/psaB etc. unflagged).

**Accession-level removal drops the real record too.** NC_034951.1 and
NC_058879.1/NC_058880.1 each carry both a correct and a mislabelled
record for the same gene (ycf3 169 aa + 198 aa; cemA 231 aa + 320 aa).
Removing the accession removes both. Dedup keeps the longest
translation, so the mislabelled one was the one selected. Follow-ups
are in tasks/todo.md.

### Integration, `v2026.10` vs `v2026.10_1` (`NXF_VER=25.04.0`, both 69/69 tasks)

- **Barcodes:** `barcodes.fasta` byte-identical for INT-ANIMAL-01,
  INT-PLANT-01-pt and INT-PLANT-01-mt. No locus lost, no coordinate
  moved.
- **INT-PLANT-01-pt:** `ycf3` 2 → 1, `cemA` 2 → 1. Other changes, all
  explained: `rpl32` 5 → 3 (two phantom copies from NC_040106.1, which
  is really rpl23, removed at 109594–109875 and 131395–131676); the
  surviving ycf3 hit is now NC_035724.1 rather than NC_031207.1 at the
  same 43913–45910 span (identity 0.9762 → 0.9765).
- **INT-ANIMAL-01:** annotation unchanged.
- **INT-PLANT-01-mt:** one `orfB` feature (NC_007886.1, identity 0.748)
  removed.
- `single_copy_genes: ["ycf3", "cemA"]` added to
  `tests/integration/expected/plant_pt/annotation_bounds.json` and
  asserted in `assertions.sh`: passes on `v2026.10_1`, fails on
  `v2026.10`; no other failures on the new output.
- Stub run passes.

### Upload

2026-10-21: user confirmed. `refdata-wf5/v2026.10_1/refs.tar.gz` and
the extracted bundle were published (`publish_refs_extracted.sh`
listing complete against the manifest). `bash scripts/fetch_refs.sh
v2026.10_1` on a clean scratch directory passed the SHA256 check and
manifest completeness check, and the fetched `proteins/` is identical
to the local build. `v2026.09_1` is left on Azure (still referenced
by branches).
