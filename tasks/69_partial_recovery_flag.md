# Task 69 — Flag partial recoveries instead of calling them "clean"

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** in progress — only the §4.3 Azure re-run remains. No preconditions.

**Related spec:** spec/03-organelles.md §3.1 (typical sizes per target),
spec/02-stages.md C6/C7 rows, spec/06a-reports.md §6a.2 and §6a.4,
spec/00-overview.md §0 (bundle layout). CONSTITUTION principle 7
(negative clarity, and "a partial recovery is not a negative"),
rules 14, 18 and 19.

## Overview

Every sample that reaches the end of the pipeline gets a
`sample_status`. For a sample that cleared the coverage gate, was
assembled, and had at least one barcode locus recovered, the status is
`ok`. `bin/collate.py` then gives it the reason **"clean recovery"**,
and the per-sample report's Outcome finding says "Clean recovery."
with a green badge.

"Clean" is too strong. The only test is "at least one locus
recovered". It does not ask how much of the organelle genome was
recovered, or how many of the panel's barcode loci were found. A real
example from 2026-10-14 (Azure run, insect sample
`barcode21_SAC1451-2-1`, `animal_mt`):

- Coverage gate `ok` at about 199×.
- The binned target was **one 5,340 bp contig**. A typical animal
  mitogenome is 15–20 kb (spec/03-organelles.md §3.1), so this is about a
  third of one.
- **COX1, COX2, COX3 and ATP6 were `not_found`.** Only CYTB and a
  105 nt ND1 passed. COX1 is the primary animal barcode.
- The status was still `ok`, "clean recovery".

An analyst skimming the run report would read that as a complete
mitogenome with a full barcode panel. CONSTITUTION principle 7 asks
for the opposite: an incomplete result must be distinguishable from a
complete one at every surface, and must still ship, because it is not
a negative.

This task adds a **recovery completeness** flag. It sits alongside
`sample_status` and does not replace it. A sample whose status is `ok`
or `low_coverage` is a **partial recovery** when either of these holds:

1. **The target is short.** The binned target's total length is under
   80% of the lower bound of the expected size range for its
   `assembly_target`.
2. **A barcode locus is missing.** Any locus in the target's panel
   (`assets/loci.json`) has status `not_found`.

Otherwise it is a **complete** recovery. `sample_status` keeps its
current six values. A partial `ok` is still `ok`, still a full bundle,
and still ships its barcodes. Only the reason text, the new metadata
block and the reports change.

## 1. Constraints

- **Principle 7.** A partial recovery is a passing result with a
  caveat. Render it as information, never as a warning or failure,
  just as `low_coverage` is rendered (spec/06a-reports.md §6a.2,
  task 54). It must also never read as complete. Both parts are
  acceptance criteria.
- **No new `sample_status` value** (decided with the user, 2026-10-14).
  The bundle dispatch (`FULL_BUNDLE_STATUSES`), `run_manifest.json`'s
  status counts, and `assertions.sh`'s `EXPECTED_SAMPLE_STATUS` stay
  as they are.
- **Rule 18 (fail closed).** If the flag cannot be computed, for
  example because the target length is unreadable or the target has no
  configured size, record that explicitly (`completeness: "unknown"`
  with a reason). Never default to `complete`.
- **Config, not code.** The sizes and the 80% fraction are pipeline
  params, set in `nextflow.config` beside `coverage_limits`, and
  validated in `main.nf`'s `validateParams()`. They are not constants
  in `collate.py`.
- **Rule 14.** All new logic in `bin/` is stdlib only and has 100%
  branch coverage.
- **Rule 19.** Integration assertions check structure and status
  (`completeness`, reason keys), not exact lengths.

## 2. Definitions

### 2.1 Expected size range

The source is spec/03-organelles.md §3.1 "Typical size". This adds a
new param, `expected_target_size`, keyed by `assembly_target`. Each
entry has `min_bp` and `max_bp`. Only `min_bp` is used by this task.
`max_bp` is recorded so the report can show the range, and so a later
"too long" check has somewhere to live.

| Target | `min_bp` | `max_bp` | 80% threshold |
|---|---|---|---|
| `animal_mt` | 15,000 | 20,000 | 12,000 bp |
| `plant_pt` | 120,000 | 160,000 | 96,000 bp |
| `plant_mt` | 200,000 | (none: "several Mb") | 160,000 bp |

`plant_mt` has no meaningful upper bound. Use `null`, and have the
report print "≥ 200 kb". Mark all three values **provisional**, in the
same way as `coverage_limits` (spec/02-stages.md §2.1.6). They come
from the literature, not from our data. Some real plant mitogenomes are
smaller than 200 kb, so a `plant_mt` flag is weaker evidence than an
`animal_mt` one.

This is deliberately separate from `coverage_limits.*.nominal_size`
(a single point estimate for the coverage arithmetic) and from
`assembly_size_hints` (a Flye hint). Don't reuse either.

New param `partial_recovery_min_fraction`, default `0.8`. It must be in
(0, 1]. The comparison is strictly less than:
`target_length_bp < fraction × min_bp` means short.

### 2.2 Target length

This is the total length, summed over all records, of the target FASTA
that COLLATE already receives (`target_fasta`, published as
`organelle_assembly.fasta`). It is **not** `assembly.total_bp`, which
is the whole metaFlye assembly. For the example sample above, that was
309 kb across 65 contigs: far larger than the mitogenome, and it would
hide the problem. For `plant_mt`, the multipartite target is the sum of
the selected contigs. For `plant_pt`, it is the canonicalised target
(spec §3.6).

### 2.3 Missing locus

A panel locus is missing if its `validation.tsv` row has status
`not_found`. These do **not** count as missing:

- A locus-level `partial` (task 51). The locus was found and a sub-span
  shipped. INT-ANIMAL-01's COX1 is one of these, and it must stay
  `complete`.
- A locus that was found but failed validation (`fail`, with a reason).
  It was found, so it is not "completely missing". The per-locus table
  already reports why it failed.

If a panel locus has no row in `validation.tsv` at all, that is an
upstream contract break, not a missing locus. Set
`completeness: "unknown"` and name the locus in the reason (rule 18).

## 3. Scope

### 3.1 `bin/collate.py`

Add a pure function beside `classify_sample()`:

```
classify_recovery(sample_status, target_length_bp, expected_min_bp,
                  min_fraction, locus_rows, panel) -> recovery dict | None

  if sample_status not in (ok, low_coverage):  -> None
      (no_barcode / no_assembly / fail / error already say what happened)
  if target_length_bp is None or expected_min_bp is None
     or a panel locus has no row:              -> completeness "unknown", reason
  short   = target_length_bp < min_fraction * expected_min_bp
  missing = [panel loci whose row status == not_found], in panel order
  completeness = "partial" if short or missing else "complete"
```

Pass the two params in as explicit CLI args (`--expected-min-bp`,
`--expected-max-bp`, `--partial-min-fraction`). The process looks them
up from `params` for `meta.assembly_target`. Don't dig them out of
`--params-json`, because explicit args are what the unit tests can
drive.

**Reason text.** Change `sample_status_reason` only when completeness
is `partial`:

- `ok` + partial: `"partial recovery — target 5,340 bp (< 80% of the
  15,000 bp expected minimum); panel loci not found: COX1, COX2, COX3,
  ATP6"`. Include only the clauses that apply.
- `low_coverage` + partial: keep the existing low-coverage reason and
  append `"; partial recovery — …"` with the same clauses.
- `ok` + complete keeps `"clean recovery"`. `ok` + unknown says
  `"recovery completeness could not be assessed — <why>"`, never
  "clean".

### 3.2 `metadata.json` and its schema

Add a top-level `recovery` object to `assets/sample_metadata.schema.json`.
It can be null (null for the statuses in §3.1's first branch). It is
not in `required`, so older bundles still validate. Example:

```
"recovery": {
  "completeness": "partial",          // complete | partial | unknown
  "reasons": ["target_short", "loci_not_found"],
  "target_length_bp": 5340,
  "expected_min_bp": 15000,
  "expected_max_bp": 20000,
  "min_fraction": 0.8,
  "loci_not_found": ["COX1", "COX2", "COX3", "ATP6"],
  "unknown_reason": null
}
```

`reasons` is a closed enum (`target_short`, `loci_not_found`), so that
the reports and the assertions key off it rather than parsing prose.
`additionalProperties: false`, as for the rest of the schema.

### 3.3 Per-sample `report.html` (`bin/report/report.py`)

- **Outcome finding:** for `ok` + partial, replace the green "Clean
  recovery." with an `info` badge. Text such as: "Partial recovery: the
  assembled target is 5,340 bp against an expected 15–20 kb, and 4 of 6
  panel loci were not found (COX1, COX2, COX3, ATP6). The recovered
  barcodes are real and usable. Read the result as incomplete, not
  failed." For `low_coverage` + partial, append the partial sentence to
  the existing low-coverage outcome. For unknown, show a `secondary`
  badge that says completeness could not be assessed.
- **Key findings:** add an "Organelle target length" row:
  `5,340 bp (expected 15–20 kb)`. It gets an `info` badge when short and
  a `success` badge otherwise. Keep it next to the existing "Assembly
  size" row, and relabel that row "Assembly size (all contigs)" so the
  two numbers can't be confused.
- Surface the partial reasons in the warnings mirror, labelled with the
  tab that owns them: Validation for missing loci, Assembly for length.

### 3.4 Run-level `run-report.html` (`bin/run_report.py`, `bin/report/run_report.py`)

- Pass `recovery` through into each row (as `sample_status_reason`
  already is) and into `run_manifest.json`'s per-sample entries
  (`assets/run_manifest.schema.json`, nullable).
- In the summary table, mark an `ok` / `low_coverage` row that is
  partial distinctly, for example an extra `partial` pill next to the
  status pill. The existing status count breakdown is unchanged. Add a
  sub-count such as `ok: 3 (1 partial)`, as spec/06a-reports.md §6a.4
  already does for locus-level partials (`5/6 (+1 partial)`).

### 3.5 Wiring

- `nextflow.config`: add `expected_target_size` and
  `partial_recovery_min_fraction`, with comments citing spec §3.1 and
  noting that the values are provisional.
- `main.nf` `validateParams()`: every assembly target has an entry,
  `min_bp` is a positive integer, `max_bp` is null or ≥ `min_bp`, and
  the fraction is in (0, 1]. Each check fails with a message naming the
  param.
- `modules/local/collate.nf`: pass the three new args.

## 4. Testing

### 4.1 Unit tests (`scripts/tests/test_collate.py`, `test_report.py`, `test_run_report.py`)

Run them with `scripts/pytest.sh`. Branch coverage for the new code
must be 100%.

- `classify_recovery`: complete; short only; missing only; both; the
  exact boundary (`length == fraction × min_bp` is **not** short); one
  under the boundary; a locus-level `partial` and a `fail` row do not
  count as missing; a panel locus absent from `validation.tsv` gives
  `unknown`; a missing target length or size config gives `unknown`; and
  `no_barcode` / `no_assembly` / `fail` give `None`.
- Reason text for `ok` × {complete, partial, unknown} and
  `low_coverage` × {complete, partial}. Check that "clean recovery"
  appears only for `ok` + complete.
- Target length summed over a multi-record FASTA (the `plant_mt` case).
- Schema: a metadata document with each `recovery` shape validates
  against `assets/sample_metadata.schema.json`. An unknown `reasons`
  value is rejected.
- Report: the Outcome finding's class and text per case. "Clean
  recovery." must never render for a partial sample. Check the target
  length row and the run-report partial marker and sub-count. Update
  `test_report.py`'s existing "Clean recovery." assertion so it still
  passes for the complete case.

### 4.2 Integration fixtures

All three current fixtures are **complete** under this definition
(checked against `tests/integration/output/` from the 2026-10-12 run):

| Fixture | Target length | Threshold | `not_found` loci |
|---|---|---|---|
| INT-ANIMAL-01 | 16,799 bp | 12,000 | none (COX1 is a locus-level `partial`) |
| INT-PLANT-01-pt | 155,276 bp | 96,000 | none |
| INT-PLANT-01-mt | 422,015 bp | 160,000 | none |

- Add an `EXPECTED_RECOVERY` map to `tests/integration/assertions.sh`
  (all three `complete`). Assert `recovery.completeness` and that
  `recovery.reasons` is empty. Also assert that INT-PLANT-01-mt's
  reason keeps the low-coverage text and has no partial clause. This
  pins the negative side and the wiring end to end.
- **No fixture exercises the partial branch**, and none is added here.
  A cheap, deterministic partial fixture (a truncated read set) would
  be fragile under assembly. Record the gap in the assertions' header
  comment, as task 42 §9 did for `no_barcode`. The partial branch is
  covered by the unit tests and by §4.3.
- Bump nothing in `expected/*/bin_bounds.json`. Those are fixture
  bounds, not the new params.

### 4.3 Real-data acceptance

Re-run `barcode21_SAC1451-2-1` on Azure with the `azure-run` skill
(`tests/manual/samples_az_barcode21_SAC1451-2-1.csv`). The work dir is
`az://workdata/wf5/20261014_063010_343359`, which is resumable until
2026-10-28, and only COLLATE onwards should re-run. Expect:

- `sample_status` `ok`, `recovery.completeness` `partial`, and reasons
  `target_short` and `loci_not_found` with COX1, COX2, COX3 and ATP6.
- The Outcome finding shows the info badge with both clauses, and the
  run report shows the partial marker.

Re-render the integration run report too. Give the user both report
paths, per CLAUDE.md.

## 5. Spec and doc updates

- spec/03-organelles.md §3.1: note that the "Typical size" lower bounds
  feed `params.expected_target_size` and the partial-recovery flag.
- spec/02-stages.md C6 row: C6 derives `recovery` (§2 of this task),
  and the "clean recovery" reason is now conditional. C7 row: carries
  `recovery` into the run report and the manifest.
- spec/06a-reports.md §6a.2 (Key findings: the target length row and
  the partial Outcome) and §6a.4 (summary-table partial marker and
  sub-count).
- spec/00-overview.md §0: one line saying that a full bundle can be
  complete or partial, and where that is recorded.
- CONSTITUTION: **no amendment needed.** Principle 7 already says a
  partial recovery is not a negative and must be distinguishable at
  every surface. This task implements that for a second axis (genome
  and panel completeness) beside coverage.
- The `azure-run` and `remote-run` skills, "Review" step: report
  `recovery.completeness` and its reasons alongside `sample_status`.

## 6. Acceptance criteria

- [x] `classify_recovery` is implemented as §2–§3.1 describe, with 100%
      branch coverage.
- [x] `metadata.json` carries `recovery`, and the schema validates it.
- [x] "Clean recovery" appears only for an `ok` sample that is complete.
- [x] Partial samples show an info (not warning or danger) Outcome in
      `report.html`, and a distinct marker in `run-report.html`.
- [x] The new params are in `nextflow.config`, validated in
      `validateParams()`, and documented as provisional.
- [x] The integration assertions pin all three fixtures as `complete`,
      and the partial-fixture gap is recorded.
- [ ] Re-running `barcode21_SAC1451-2-1` gives `partial` with both
      reasons. Report paths are given to the user.
- [x] The spec and skill updates in §5 are made.


## Outcomes

- Implemented as specified: `classify_recovery`, `apply_recovery_reason`,
  `fasta_total_length` and `panel_loci` in `bin/collate.py`; `recovery`
  in both schemas; target-length row, Outcome variants and warnings in
  `bin/report/report.py`; `partial` pill and sub-count in the run
  report; params in `nextflow.config`, validated in `validateParams()`.
  Full suite: 961 passed, 100% branch coverage on `bin/`; flake8 clean.
- **Deviations / additions (small, needed to make the spec work):**
  - COLLATE also takes `--locus-panel` (`assets/loci.json`, staged from
    `params.locus_panel`), since the panel is required to detect a locus
    with no `validation.tsv` row. `classify_recovery` also takes
    `expected_max_bp`, so it can sit in the `recovery` block.
  - A locus with several `validation.tsv` rows (e.g. plant_mt `nad1`) is
    missing only if every row is `not_found`.
  - An empty/absent panel gives `unknown`, not `complete` (rule 18).
  - The Outcome sentence says "N panel loci were not found" without
    "of M", as `recovery` does not carry the panel size.
- **Check against the 2026-10-14 outputs (local, no re-run):**
  `classify_recovery` on `tests/manual/output/az_barcode21_SAC1451-2-1/`
  gives `partial`, reasons `target_short` + `loci_not_found`, loci COX1,
  COX2, COX3, ATP6, target 5,340 bp, and the expected reason text.
- **Not done: §4.3 Azure re-run.** `tests/manual/samples_az_barcode21_SAC1451-2-1.csv`
  does not exist and the reads' `az://` URI is not in the local logs.
  The two report paths for the user are pending that run.
