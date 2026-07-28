# Task 4 — P1 stage 0: `PARSE_SAMPLESHEET`

> **Status:** Shipped against the initial `kingdom` column. The samplesheet
> schema has since been redesigned to use `assembly_target` (one row → one
> organelle — see [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)).
> Migration is tracked in [task 9](9_assembly_target_migration.md). The
> validation cases below still apply verbatim, substituting
> `assembly_target ∈ {animal_mt, plant_pt, plant_mt}` wherever this task
> says "bad kingdom" / `kingdom`.

**Phase:** P1 (first real workflow stage — from [plan.md §6](../plan.md)).
**Goal:** Replace the P0 stub for stage 0 with the real samplesheet-parsing
logic. This is [C1 in spec §2.2](../spec/02-stages.md#22-custom-logic-components):
CSV validation beyond nf-schema's reach, plus per-sample read concatenation.

**Exit criteria:**

- `nextflow run main.nf -profile test` produces, for every row in
  `tests/data/samples.csv`, a single concatenated
  `results/<sample_id>/qc/<sample_id>.reads.fastq.gz` (when
  `publish_intermediates=true`) whose contents are the byte-concatenation
  of the resolved input FASTQs in samplesheet order.
- Invalid samplesheets fail preflight with a clear, actionable error and
  a non-zero exit **before** any per-sample process is scheduled. The
  error names the offending row and the offending field. Covered cases
  (see §3 below): duplicate `sample_id`, absolute path in `reads`,
  missing FASTQ file, bad `kingdom` value, unknown column in header.
- Unit tests under `scripts/tests/` cover each validation branch and pass
  under `pytest` inside the `neoformit/daff-wf5-scripts` container.
- No `${params.<file>}` bare interpolation regression
  ([plan.md §1a](../plan.md)); all file staging goes through explicit
  channels or `${file(params.x)}`.

**Not in scope:**

- Any change to QC / recruitment / downstream stages — they still consume
  `PARSE_SAMPLESHEET.out.reads` with the same tuple shape as P0.
- Uploading a rebuilt `neoformit/daff-wf5-scripts` image — that ships via
  the CI `images.yml` trigger once `requirements.txt` or the Dockerfile
  changes ([spec §2.2](../spec/02-stages.md#22-custom-logic-components),
  [spec §5b](../spec/05-test-data.md)).

**Cross-cutting rules (from [plan.md §1a](../plan.md)):**

- Runs in `neoformit/daff-wf5-scripts:<tag>` (the shared bespoke image
  built from [`scripts/Dockerfile`](../scripts/Dockerfile) +
  [`scripts/requirements.txt`](../scripts/requirements.txt)) — no host
  Python, no conda.
- Python source is **not** baked into the image. It lives under `bin/`
  and is auto-staged onto the container `PATH` by Nextflow at runtime
  ([spec §2.2](../spec/02-stages.md#22-custom-logic-components) "deps in
  the image, code at runtime").
- Every file crossing a process boundary is a Nextflow `Path`. The
  samplesheet is materialised once as
  `ch_samplesheet = Channel.value(file(params.samplesheet))`; resolved
  FASTQs enter per-sample channels as `path` inputs.

---

## 1. Architecture

Two-part split — the same one flagged in
[spec §2.2 C1](../spec/02-stages.md#22-custom-logic-components):

1. **`VALIDATE_SAMPLESHEET`** — a single preflight process that runs
   *once per invocation*, before any fan-out. It validates the CSV
   holistically (uniqueness, existence, kingdom enum, header shape) and
   emits a normalised JSON array — one object per validated sample, with
   reads resolved to absolute paths. Fails the whole run with a non-zero exit
   and a clear error message if anything is wrong. This replaces the
   Groovy `buildMeta` / `resolveReads` helpers currently in
   [`main.nf`](../main.nf) — validation belongs in Python + tests, not in
   an untested `.nf` closure.
2. **`PARSE_SAMPLESHEET`** — per-sample process (already stubbed in
   [`modules/local/parse_samplesheet.nf`](../modules/local/parse_samplesheet.nf)).
   Takes `tuple(meta, path(reads_list))` and produces the concatenated
   `${meta.sample_id}.reads.fastq.gz` that feeds `NANOPLOT_RAW` and
   `CHOPPER`. Single-file rows are a passthrough (still re-emitted under
   the canonical filename for downstream tag consistency); multi-file
   rows are `cat`-concatenated in samplesheet order.

Rationale for the split: preflight failure must abort the whole run
before scheduling anything (a bad samplesheet is a user error, not a
per-sample soft-fail). Per-sample concat only needs the resolved reads
list, so it fans out cleanly.

### 1.1 `main.nf` shape after this task

```groovy
ch_samplesheet = Channel.value(file(params.samplesheet))
ch_data_dir    = Channel.value(file(params.data_dir))

VALIDATE_SAMPLESHEET(ch_samplesheet, ch_data_dir)

ch_samples = VALIDATE_SAMPLESHEET.out.json
    | splitJson()
    | map { row ->
        def meta = [
            sample_id          : row.sample_id,
            kingdom            : row.kingdom,
            sample_info        : row.sample_info         ?: '',
            sample_type        : row.sample_type         ?: '',
            sample_receipt_date: row.sample_receipt_date ?: '',
            storage_location   : row.storage_location    ?: '',
        ]
        // reads are already absolute, resolved, existence-checked
        def reads = row.reads.collect { file(it) }
        [ meta, reads ]
    }

PARSE_SAMPLESHEET(ch_samples)
```

Delete `validateParams()`, `buildMeta()`, `resolveReads()` from
[`main.nf`](../main.nf) — their responsibilities move into the Python
validator, which is testable in isolation.

## 2. `bin/parse_samplesheet.py`

Single entry point invoked from `VALIDATE_SAMPLESHEET`. CLI:

```
parse_samplesheet.py --samplesheet <csv> --data-dir <dir> --out <json>
```

Exits 0 on success (JSONL written); non-zero with a human-readable
error on any validation failure. All errors printed to stderr,
prefixed with `ERROR:` and naming `<row N> (sample_id=<id>)`. On success,
writes a single JSON array to `--out` (one object per row, keys in the
order below).

Behaviour, per [spec §0](../spec/00-overview.md) validation rules:

| Check | Failure mode |
|---|---|
| Header has `sample_id`, `kingdom`, `reads` | hard fail — missing required column |
| Header contains only known columns (required ∪ optional) | hard fail — unknown column `<name>` (guards against silent typo drops) |
| `sample_id` matches `^[A-Za-z0-9_.-]+$` | hard fail — bad sample_id |
| `sample_id` unique across the sheet | hard fail — duplicate sample_id `<id>` |
| `kingdom` ∈ {`plant`, `animal`} (case-insensitive → lowercased) | hard fail — bad kingdom |
| `reads` split on `\|`; each piece non-empty | hard fail — empty reads entry |
| No `reads` path starts with `/` | hard fail — "absolute path not allowed; use --data-dir root" |
| Each resolved read path exists under `--data-dir` | hard fail — file not found (prints full resolved path) |
| Each resolved path ends in `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz` | hard fail — unsupported extension |
| `sample_receipt_date`, if present, parses as ISO 8601 `YYYY-MM-DD` | **warn only** — pass raw string through; log to stderr |
| Optional columns absent from header | ok — treat as empty for every row |

Output JSON — a single array of objects, all optional keys always present:

```json
[
  {
    "sample_id": "TEST-PLANT-01",
    "kingdom": "plant",
    "reads": ["/abs/path/to/plant.fastq.gz"],
    "sample_info": "synthetic plant test fixture",
    "sample_type": "leaf",
    "sample_receipt_date": "2026-07-24",
    "storage_location": "N/A"
  }
]
```

Absolute paths in the JSON let the downstream `.map { file(it) }` in
`main.nf` produce stageable `Path` objects on any executor
([plan.md §1a](../plan.md)).

Implementation notes:

- stdlib `csv` + `json` are enough; no need to pull `pandas` into the
  hot path for a file that fits in memory. `pyyaml` is already in
  `requirements.txt` for other consumers, so no new deps.
- `datetime.date.fromisoformat` for the date warn-only check.
- Collect **all** validation errors before exiting — a user with 30
  malformed rows should see all 30, not fix-and-retry 30 times. Exit
  once at the end with the accumulated message.

## 3. `modules/local/validate_samplesheet.nf` (new)

```groovy
process VALIDATE_SAMPLESHEET {
    label 'process_low'

    input:
    path samplesheet
    path data_dir

    output:
    path 'samples.normalised.json', emit: json

    script:
    """
    parse_samplesheet.py \\
        --samplesheet ${samplesheet} \\
        --data-dir ${data_dir} \\
        --out samples.normalised.json
    """
}
```

No `tag`; no per-sample fan-out at this point. `container` resolved via
`conf/containers.config` — add a `VALIDATE_SAMPLESHEET` entry pointing
at the same `neoformit/daff-wf5-scripts:<tag>` used by
`PARSE_SAMPLESHEET`.

## 4. `modules/local/parse_samplesheet.nf` (replace stub)

```groovy
process PARSE_SAMPLESHEET {
    tag        "${meta.sample_id}"
    label      'process_low'
    publishDir "${params.outdir}/${meta.sample_id}/qc", mode: 'copy',
               enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.sample_id}.reads.fastq.gz"), emit: reads

    script:
    """
    # Byte-concat of gzipped FASTQs is a valid gzipped stream; no re-gzip.
    # Single-file rows still go through cat so the output filename is canonical.
    cat input/* > ${meta.sample_id}.reads.fastq.gz
    """
}
```

Notes:

- `stageAs: 'input/*'` avoids filename collisions when a row references
  the same basename twice from different subdirs (rare but not
  impossible; better to be defensive).
- `cat` on concatenated `.gz` streams is standard — `zcat` / gzip
  readers reassemble them transparently. No need to decompress and
  re-compress here.
- Container: reuse `neoformit/daff-wf5-scripts:<tag>` (has `coreutils`),
  or a `coreutils` biocontainer if we want to keep the scripts image
  lean. Recommend reusing the scripts image — one fewer image to pull.

## 5. `conf/containers.config`

Add / update:

```groovy
withName: 'VALIDATE_SAMPLESHEET' {
    container = 'neoformit/daff-wf5-scripts:<sha>'
}
withName: 'PARSE_SAMPLESHEET' {
    container = 'neoformit/daff-wf5-scripts:<sha>'
}
```

Replace the P0 `python:3.12-slim` placeholders. `<sha>` = the first
image tag built from the current `scripts/Dockerfile` +
`requirements.txt` — no `latest` ([plan.md §1a](../plan.md)).

## 6. Tests

### 6.1 Unit — `scripts/tests/test_parse_samplesheet.py`

Under `pytest`, one test class per validation branch. Fixtures live in
`scripts/tests/data/parse_samplesheet/` — each is a small CSV + a
minimal fake `--data-dir` tree of empty FASTQs.

Cases (each maps to a table row in §2):

- `test_valid_multirow_sheet_emits_jsonl_in_order`
- `test_valid_multifile_row_resolves_all_reads`
- `test_duplicate_sample_id_fails`
- `test_absolute_reads_path_fails`
- `test_missing_reads_file_fails`
- `test_bad_kingdom_fails`
- `test_kingdom_case_normalised_to_lowercase`
- `test_unknown_column_in_header_fails`
- `test_missing_required_column_fails`
- `test_bad_sample_id_pattern_fails`
- `test_unsupported_reads_extension_fails`
- `test_bad_date_warns_but_passes_through`
- `test_all_errors_collected_before_exit` — assert that a sheet with 3
  distinct problems reports all 3 in the stderr message.
- `test_optional_columns_absent_from_header_ok`
- `test_optional_columns_empty_cells_ok`

Run inside the `neoformit/daff-wf5-scripts:<sha>-test` container in CI
(matches runtime environment plus `pytest`). Add `pytest` to
`scripts/requirements-test.txt`; `scripts/Dockerfile` picks it up when
built with `--build-arg TEST=1` (see §8).

### 6.2 Integration — `-profile test`

Extend `tests/data/samples.csv` to include one multi-file row so the
concat path is exercised end-to-end:

```
sample_id,kingdom,reads,sample_info,sample_type,sample_receipt_date,storage_location
TEST-PLANT-01,plant,plant.fastq.gz,synthetic plant test fixture,leaf,2026-07-24,N/A
TEST-ANIMAL-01,animal,animal.fastq.gz,synthetic animal test fixture,whole,2026-07-24,N/A
TEST-ANIMAL-02,animal,animal.fastq.gz|animal.fastq.gz,two-file concat test,whole,2026-07-25,N/A
```

Assert (in a shell test step after `nextflow run`): the concat output
for `TEST-ANIMAL-02` is exactly `2 × size(animal.fastq.gz)` bytes.

Add a negative-path integration test — a broken samplesheet fixture at
`tests/data/samples_bad.csv` (duplicate `sample_id`), run under a
`test_bad_samplesheet` profile that points at it, and assert `nextflow
run` exits non-zero **before** any per-sample process runs (no
`work/<hash>` dirs for `NANOPLOT_RAW` etc.). This proves the preflight
gate works.

## 7. Deliverables checklist

- [ ] `bin/parse_samplesheet.py` implemented, `flake8`-clean.
- [ ] `modules/local/validate_samplesheet.nf` added.
- [ ] `modules/local/parse_samplesheet.nf` replaced with the real
      concat implementation.
- [ ] `main.nf` refactored: `buildMeta`/`resolveReads`/`validateParams`
      removed; `VALIDATE_SAMPLESHEET` wired in as preflight;
      `splitJson()` feeds the per-sample channel.
- [ ] `conf/containers.config` — `VALIDATE_SAMPLESHEET` +
      `PARSE_SAMPLESHEET` pointed at
      `neoformit/daff-wf5-scripts:<sha>` (pinned tag).
- [ ] `scripts/tests/test_parse_samplesheet.py` covers every branch in
      §6.1; passes under the scripts container.
- [ ] `tests/data/samples.csv` extended with the multi-file row;
      `tests/data/samples_bad.csv` added.
- [ ] `-profile test` passes with the extended fixture; `-profile
      test_bad_samplesheet` fails preflight cleanly.
- [ ] CI green — lint + container-coverage + `-profile test` +
      `pytest` on `scripts/tests/`.

## 8. Resolved design decisions

1. **Keep nf-schema alongside the Python validator.** nf-schema handles
   header/type shape (fast fail, no container spin-up); the Python
   validator handles cross-row + filesystem checks nf-schema cannot
   express (uniqueness, existence, absolute-path rejection). The
   overlap is small and the layered defence is cheap.
2. **Test deps live in a separate `scripts/requirements-test.txt`.**
   The runtime image stays lean; CI builds a `-test` variant that
   layers `pytest` on top. `scripts/Dockerfile` gains an optional
   `TEST=1` build arg that installs the extra file.
