# Task 9 — Migrate samplesheet from `kingdom` to `assembly_target`

**Phase:** P1 refactor. Prerequisite for [task 8](8_recruit.md) (RECRUIT
real implementation) and every downstream stage whose parameters key on
the target organelle.
**Goal:** Migrate the codebase from the shipped `kingdom` column
(plant / animal) to the redesigned `assembly_target` column
(`animal_mt` / `plant_pt` / `plant_mt`) per the cross-cutting rule in
[spec §1a](../../spec/01-pipeline-flow.md#1a-engineering-constraints).
The design docs already reflect the new schema; this task lands the
code changes so they stop drifting.

**Exit criteria:**

- `nextflow run main.nf -profile test` passes 47/47 against a
  redesigned `tests/data/samples.csv` where the `kingdom` column is
  replaced with `assembly_target`, using the fixed enum
  `{animal_mt, plant_pt, plant_mt}`.
- `nextflow run main.nf -profile test_bad_samplesheet` still exits
  non-zero at preflight, and the negative-path fixture still triggers
  a validation error.
- `meta.assembly_target` is threaded through every process input tuple.
  `meta.kingdom` no longer exists.
- All 17 pytest cases in
  [`scripts/tests/test_parse_samplesheet.py`](../../scripts/tests/test_parse_samplesheet.py)
  pass after being updated to reference `assembly_target`. New test
  case added: rejecting a row with an unknown `assembly_target` value
  (e.g. `plant`, `animal`, `fungal_mt`).
- CI green — lint + container-coverage + `-profile test` + `pytest`.

**Not in scope:**

- Any real-tool wiring for stages downstream of PARSE_SAMPLESHEET —
  those land in their own tasks ([task 8](8_recruit.md) onward). This
  task is a pure column rename + enum widening; behaviour of existing
  QC stages (CHOPPER, FILTLONG, NanoPlot) is unaffected because they
  don't consume the `kingdom` / `assembly_target` field.
- Renaming `params.kingdom_refs` → `params.organelle_refs` is **in
  scope** even though it's the *reference* param, not the samplesheet
  column, because leaving them out of sync creates lasting confusion.
- Rebuilding the reference bundle at `params.organelle_refs`. That
  bundle already lays out per-target `.mmi` files
  ([task 3 §1](3_refdata.md#1-target-layout)) — no bundle change is
  required.

**Cross-cutting rules (from [spec §1a](../../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- One sample-sheet row → one organelle. Every stage keys off
  `meta.assembly_target`.
- `meta.kingdom` (if still needed anywhere for reporting) is derived
  from `assembly_target` (`plant` if it starts with `plant_`, else
  `animal`), never read from the sheet.

---

## 1. Files to change

### 1.1 Python validator

[`bin/parse_samplesheet.py`](../../bin/parse_samplesheet.py)

- Rename constant `VALID_KINGDOMS = {'plant', 'animal'}` →
  `VALID_ASSEMBLY_TARGETS = {'animal_mt', 'plant_pt', 'plant_mt'}`.
- In `_validate_row`, rename local `kingdom_raw` / `kingdom` →
  `assembly_target_raw` / `assembly_target`. The `.lower()`
  normalisation and enum check stay identical in shape.
- Emit `"assembly_target"` (not `"kingdom"`) in the output JSON object.
- Update the row-labelling error string to name
  `assembly_target` in the "invalid" message.

### 1.2 Nextflow workflow

[`main.nf`](../../main.nf)

- In the `splitJson()` → `map` closure, build `meta` with
  `assembly_target: row.assembly_target`; drop `kingdom`. Do not
  synthesise a `meta.kingdom` field — downstream that needs it can
  derive from the target (`meta.assembly_target.startsWith('plant_')`).
- Rename `params.kingdom_refs` → `params.organelle_refs` throughout
  main.nf (channel variable `ch_kingdom_refs` → `ch_organelle_refs`).

### 1.3 Nextflow config

[`nextflow.config`](../../nextflow.config)

- Rename `params.kingdom_refs` → `params.organelle_refs` in the
  `params { … }` block.

`conf/test.config`,
`conf/test_bad_samplesheet.config`

- Same rename. Value continues to point at the placeholder
  `tests/data/refs/kingdom_refs.mmi` (or the empty file) until
  [task 8](8_recruit.md) lands and swaps in the recruit directory.

### 1.4 Process modules

Every module that takes an input named `kingdom_refs` renames it to
`organelle_refs`. As of writing, that's:
- [`modules/local/recruit.nf`](../../modules/local/recruit.nf) (stub)
- [`modules/local/bin_target.nf`](../../modules/local/bin_target.nf) (stub)

No script-block behaviour change — these are stubs. The rename is
purely cosmetic for now, but locks in the naming convention downstream
tasks will use.

### 1.5 Schema

[`assets/samplesheet.schema.json`](../../assets/samplesheet.schema.json)
was already updated in the design pass — no further change needed.
Verify the schema matches `{animal_mt, plant_pt, plant_mt}` and
`assembly_target` is in `required`.

### 1.6 Fixtures

`tests/data/samples.csv`

Redesigned rows. Note that `TEST-PLANT-01` was previously a single
plant row; under the new schema it becomes either a single
`plant_pt` row or two rows if we want to test both organelles.

Recommend two plant rows to exercise the multi-row-same-reads pattern
end-to-end:

```csv
sample_id,assembly_target,reads,sample_info,sample_type,sample_receipt_date,storage_location
TEST-PLANT-01-pt,plant_pt,plant.fastq.gz,synthetic plant test fixture,leaf,2026-07-24,N/A
TEST-PLANT-01-mt,plant_mt,plant.fastq.gz,synthetic plant test fixture,leaf,2026-07-24,N/A
TEST-ANIMAL-01,animal_mt,animal.fastq.gz,synthetic animal test fixture,whole,2026-07-24,N/A
TEST-ANIMAL-02,animal_mt,animal.fastq.gz|animal_b.fastq.gz,two-file concat test,whole,2026-07-25,N/A
```

Update `tests/data/samples_bad.csv`
similarly — keep the duplicate `sample_id` failure case, and swap
`kingdom` column for `assembly_target`.

### 1.7 Unit tests

[`scripts/tests/test_parse_samplesheet.py`](../../scripts/tests/test_parse_samplesheet.py)

- Global find/replace `kingdom` → `assembly_target` in test CSV fixtures.
- Update `test_bad_kingdom_fails` → `test_bad_assembly_target_fails`;
  replace fixture value `fungi` with `unknown_target` (or
  `plant`, since a bare kingdom without organelle is now invalid).
- Update `test_kingdom_case_normalised_to_lowercase` → same for
  `assembly_target` (e.g. `Plant_Pt` → `plant_pt`).
- Add `test_kingdom_column_rejected` — a CSV with a legacy `kingdom`
  column (and no `assembly_target`) is rejected with the standard
  "unknown column" or "missing required column" error. Guards against
  silent regression when a stale samplesheet is submitted.
- Assert output JSON contains `assembly_target`, not `kingdom`.

### 1.8 Documentation

`tests/data/README.md` — update the sample
sheet section to describe the new column and the multi-row pattern.

## 2. Migration order

Sequence matters — do these in one commit to avoid a broken
intermediate state:

1. Python validator + unit tests (verify locally with `pytest`).
2. Schema JSON (already done in design pass — verify still correct).
3. `main.nf` meta building + param rename.
4. `nextflow.config` param rename.
5. `conf/test.config` + `conf/test_bad_samplesheet.config` param
   rename.
6. Process modules `kingdom_refs` → `organelle_refs` input rename.
7. `tests/data/samples.csv` + `samples_bad.csv` fixture update.
8. `tests/data/README.md` doc update.
9. Run `-profile test` and `-profile test_bad_samplesheet`; confirm
   both behave as expected.

## 3. Verification

- `pytest scripts/tests/test_parse_samplesheet.py -v` — all cases pass,
  including new legacy-column-rejection case.
- `nextflow run main.nf -profile test` — 47/47 pass; per-sample
  directories under `results/` still populated.
- `nextflow run main.nf -profile test_bad_samplesheet` — exits non-zero
  at VALIDATE_SAMPLESHEET (duplicate sample_id still fires).
- `grep -rn '\bkingdom\b' main.nf bin/ modules/ conf/ scripts/tests/` —
  returns only intentional derived-value references (e.g. a comment
  noting `plant_*` implies plant kingdom). Zero live column reads.
- `grep -rn 'kingdom_refs' main.nf conf/ modules/ nextflow.config` —
  returns nothing.

## 4. Deliverables checklist

- [ ] `bin/parse_samplesheet.py` — `assembly_target` throughout, new
      enum, output JSON key renamed.
- [ ] `main.nf` — meta uses `assembly_target`; `params.organelle_refs`
      + `ch_organelle_refs` rename.
- [ ] `nextflow.config` — `params.organelle_refs`.
- [ ] `conf/test.config`, `conf/test_bad_samplesheet.config` — same
      rename.
- [ ] `modules/local/{recruit,bin_target}.nf` — input rename
      `kingdom_refs` → `organelle_refs`.
- [ ] `tests/data/samples.csv` — `assembly_target` column, plant rows
      split into `plant_pt` + `plant_mt`.
- [ ] `tests/data/samples_bad.csv` — column renamed.
- [ ] `scripts/tests/test_parse_samplesheet.py` — 17 existing cases
      updated + one new legacy-column-rejection case.
- [ ] `tests/data/README.md` — column doc updated.
- [ ] `tasks/completed.txt` — task 9 appended.
- [ ] CI green.

## 5. Notes / non-issues

- **Downstream tasks (task 8, task 10+) assume `meta.assembly_target`
  exists.** They already reference it in their own specs. This task is
  the enabling change.
- **Legacy samplesheets** — no backward compatibility. A user
  submitting a sheet with `kingdom` gets a clear "missing required
  column: assembly_target" error. Documented in the migration note at
  the top of [task 4](4_parse_samplesheet.md).
- **`meta.kingdom` derivation** — deliberately not synthesised into
  `meta` at parse time. Any downstream stage that needs kingdom
  (e.g. for a report grouping) derives inline:
  `def kingdom = meta.assembly_target.startsWith('plant_') ? 'plant' : 'animal'`.
  Keeps `meta` clean and avoids a shadow-of-truth field.
- **Reference bundle filenames must match the enum values exactly.**
  Task 3 layout uses `plant_pt.mmi` / `plant_mt.mmi` / `animal_mt.mmi`
  — 1:1 with the schema enum. This task doesn't touch the bundle;
  [task 8](8_recruit.md) verifies the match at RECRUIT time.
