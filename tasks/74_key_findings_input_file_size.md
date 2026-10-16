# Task 74 — "Input file size" row in Key findings

## Overview

The Key findings table at the top of the per-sample report's Overview tab
tells the reader what the pipeline produced: contig count, assembly size,
estimated coverage, annotated genes, barcodes recovered, the top BLAST hit
and the overall outcome. It says nothing about how much data went in.

The size of the upload is the first thing an analyst checks when a result
looks thin. A 4 MB upload that fails the coverage gate is a sequencing
problem. A 3 GB upload that fails it suggests the wrong kingdom was
declared, or that the target organelle is a tiny share of the reads.
Today the reader has to find the uploaded files and check their size by
hand.

This task adds an **"Input file size"** row to Key findings. It shows the
total on-disk size of the read files submitted for the sample, in a
human-readable form such as `112 MB` or `1.45 GB`. No pipeline stage
changes: the per-file byte counts are already measured at intake (task
61) and already inlined into `metadata.json`. This is a report-only
change.

## 1. Data source

`metadata.json` → `intake` (the verbatim `intake.json` from
PARSE_SAMPLESHEET, schema `wf5/intake/v1`). Each entry in
`intake.files[]` carries `bytes`, the file's on-disk size
(`path.stat().st_size` in `bin/intake_reads.py`), so it is the
compressed size for `.fastq.gz` and the plain size for `.fastq`. Example
from `tests/integration/output/INT-ANIMAL-01/qc/intake.json`:

```
"files": [
  {"name": "animal_mt.fastq.gz", "compressed": true,
   "records": 1242, "bytes": 7541164, "error": null}
]
```

The row's value is the **sum of `bytes` over `intake.files[]`**. A
sample's reads can come from several files (task 56, a directory or a
list of read sources), and the reader cares about the total upload.

Files that failed the integrity check still have a `bytes` value (the
stat happens whatever the scan finds), so the figure is there for
`invalid_input` samples too.

No schema change is needed. `intake.files` is already in
`assets/sample_metadata.schema.json` as an untyped array.

## 2. Formatting

Add a small helper to the report package (e.g. `bin/report/utils.py`)
that turns a byte count into a human-readable string. Keep it inside
`bin/report/`. Do not import `human_bytes` from `bin/intake_reads.py`:
that function only covers MB/GB and is meant for error messages. Write
the units as module constants, following the house style.

- Decimal (SI) units: 1 KB = 1,000 bytes, 1 MB = 10^6, 1 GB = 10^9. This
  matches `bin/intake_reads.py` and the sizes the facility's storage
  shows.
- Three significant figures, with a space before the unit. The space
  matches how the report already writes other units (`15,340 bp`).

| bytes | display |
|---|---|
| 512 | `512 B` |
| 7,541,164 | `7.54 MB` |
| 112,400,000 | `112 MB` |
| 1,452,000,000 | `1.45 GB` |
| 12,300,000,000 | `12.3 GB` |

Pseudocode:

```
for unit, factor in [(GB, 1e9), (MB, 1e6), (KB, 1e3)]:
    if n >= factor: return f"{3-significant-figure value of n/factor} {unit}"
return f"{n} B"
```

Watch the rounding edge: 999,950,000 bytes must not render as
`1000 MB`. Promote it to `1.00 GB`, or pick the unit after rounding.

**File count.** When `intake.files` holds more than one file, add the
count: `1.45 GB (12 files)`. With one file, show only the size.

## 3. The row

In `bin/report/report.py`, add a `_input_size_finding(metadata)` helper
next to the other `_…_finding` helpers. It returns:

- `label`: `Input file size`
- `text`: the formatted total from §2, or `-` when `intake` is null,
  has no `files`, or no file carries an integer `bytes`. Older bundles
  and runs without `--intake-json` have no intake block, and they must
  render the same "not available" dash the other rows use, never `0 B`.
- `class`: `secondary`, the neutral badge "Assembly size" uses. Input
  size is context, not a verdict. The coverage gate has already turned
  the data volume into a pass, warn or fail outcome, and a size
  threshold here would compete with it.

**Placement.** Input size is an input property, so it comes before the
result rows. Negative clarity still wins (CONSTITUTION.md principle 7;
spec/06a-reports.md §6a.1): on a terminal page the
danger Outcome row stays first.

- `ok` / `low_coverage` / `no_barcode`: Input file size is the **first**
  row, before "Assembled contigs".
- `invalid_input` / `fail` / `no_assembly`: Outcome row first, Input file
  size second. On these pages the size is often the quickest explanation
  of what went wrong, e.g. a tiny upload behind a coverage-gate fail. It
  is worth showing even though no assembly rows follow.

The template needs no change if it already iterates the findings list
generically. Check `bin/report/templates/` and confirm the new row picks
up the standard badge-before-value layout (task 54 §3), not the compact
font used for the top BLAST hit.

## 4. Spec update

In spec/06a-reports.md §6a.2, tab 1, block **(b) Key findings**, add one
sentence: the table opens with an "Input file size" row (total on-disk
size of the submitted read files, from `metadata.json` `intake`, neutral
badge, file count when more than one). On terminal pages it follows the
Outcome row. Cite "task 74_key_findings_input_file_size.md" in plain
text.

## 5. Testing

### Unit tests (`scripts/tests/test_report.py`, run via `scripts/pytest.sh`)

Formatter:
- One case per row of the §2 table.
- 0 bytes → `0 B`.
- The 999,950,000 rounding edge renders in GB, not as `1000 MB`.

`_input_size_finding` / `key_findings`:
- One file → size only, no file count.
- Several files → summed bytes and `(N files)`.
- `intake` is null, `files` is empty, or `bytes` is missing → text `-`,
  class `secondary`.
- `ok`, `low_coverage`, `no_barcode` → `findings[0]["label"] ==
  "Input file size"`.
- `invalid_input`, `fail`, `no_assembly` → `findings[0]` is still the
  danger Outcome row and `findings[1]` is Input file size.
- Update `test_invalid_input_key_finding_is_distinct`, which asserts
  `len(findings) == 1`. It should now expect two rows and still check
  that the first is the danger Outcome naming the corrupt file.
- Update `base_metadata()` / `_sparse()` fixtures to include an
  `intake` block wherever existing tests index findings by position, so
  they keep passing for the right reason.
- A render test: `Input file size` and the formatted value appear inside
  `id="key-findings"`.

Keep coverage at or above the current gate (task 15).

### Integration

The workflow is unchanged, so no fixtures need regenerating:
`metadata.json` already carries `intake` in every integration output.
Add one assertion to `tests/integration/assertions.sh` that each
non-error sample's `report.html` contains `Input file size`. That
catches a broken wiring between `metadata.json` and the report.

### Lint

`/home/cameron/.local/envs/claude/bin/flake8 bin/report/ scripts/tests/`

## 6. Done when

- The Key findings table on every per-sample report shows "Input file
  size", placed as §3 describes, with a correctly formatted value or `-`.
- Unit tests above pass. Coverage gate holds. flake8 is clean.
- spec/06a-reports.md is updated as §4 describes.
- A report has been re-rendered (e.g. from
  `tests/integration/output/INT-ANIMAL-01/metadata.json`) and its path
  given to the user for visual inspection.

## References

- spec/06a-reports.md §6a.1, §6a.2 (Overview → Key findings), §6a.2a
- task 61_intake_integrity_and_failure_visibility.md (source of
  `intake.json` and per-file `bytes`)
- task 56_samplesheet_read_sources.md (multiple read files per sample)
- task 54_report_ui_enhancements.md §3 (Key findings badge layout)
- task 69_partial_recovery_flag.md also adds a Key findings row
  ("Organelle target length"). If both land, the ordering is Input file
  size → Assembled contigs → Assembly size (all contigs) → Organelle
  target length → …
