# Development Plan: Organelle Genome Assembly + Barcode Recovery

**Status:** Draft v0.2 · **Companion:** [brief.md](brief.md) · **Format:** Nextflow DSL2

This document expands [brief.md](brief.md) §6 into a concrete stage-by-stage
pipeline plan and lays out the development phases. Where this plan and the
brief disagree, this plan reflects the more recent decisions and the brief
should be reconciled (see §6).

---

## 0. Input: sample sheet

Pipeline is invoked with **two required parameters**:

- `--samplesheet <path>` — path to `samples.csv` (see schema below).
- `--data-dir <path>` — root directory that all `reads` paths in the sheet are resolved against.

One row per sample; multiple samples per invocation are supported and run in
parallel. All downstream stages fan out per-sample; outputs are keyed by
`sample_id`.

**Schema:**

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Unique per row. Used as the directory name for per-sample outputs. Must match `[A-Za-z0-9_.-]+` — no whitespace, no path separators. |
| `kingdom` | yes | `plant` \| `animal`. Rows with any other value (including empty) are rejected at parse time ([brief.md §1](brief.md)). |
| `reads` | yes | One or more **relative paths** to ONT FASTQ files, **pipe-delimited**. Paths are resolved against `--data-dir` (e.g. sheet says `run17/A.fastq.gz\|run17/A_pass2.fastq.gz`, `--data-dir /data`, → `/data/run17/A.fastq.gz` etc.). Absolute paths in the sheet are rejected. Multiple files are concatenated for that sample before `NANOPLOT_RAW`. `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz` all accepted. |
| `sample_info` | no | Free-text sample information. Carried through to per-sample `metadata.json` and rendered in `report.html`. Ignored by pipeline logic. |
| `sample_type` | no | Sample type descriptor (submitter-defined vocabulary). Carried through to `metadata.json` and `report.html`. Ignored by pipeline logic. |
| `sample_receipt_date` | no | Date the sample was received. Accepted format: ISO 8601 `YYYY-MM-DD`. Parse-warned (not failed) on other formats; the raw string is passed through. |
| `storage_location` | no | Physical storage location of the sample (submitter-defined). Carried through to `metadata.json` and `report.html`. Ignored by pipeline logic. |

**Example:**

Invocation: `nextflow run main.nf --samplesheet samples.csv --data-dir /data/run17`

```csv
sample_id,kingdom,reads,sample_info,sample_type,sample_receipt_date,storage_location
INT-2026-0007,plant,INT-2026-0007.fastq.gz,dried leaf fragment,leaf,2026-07-14,Freezer B / shelf 3 / box 12
INT-2026-0008,animal,INT-2026-0008_a.fastq.gz|INT-2026-0008_b.fastq.gz,merged two flowcells,,whole specimen,2026-07-18,
```

Resolved: `/data/run17/INT-2026-0007.fastq.gz`,
`/data/run17/INT-2026-0008_a.fastq.gz`, `/data/run17/INT-2026-0008_b.fastq.gz`.

**Validation (stage 0):**

- Header row: required columns (`sample_id`, `kingdom`, `reads`) must be present; optional columns (`sample_info`, `sample_type`, `sample_receipt_date`, `storage_location`) are recognised if present but not required. Order-flexible; case-sensitive. Unknown columns fail parse (guards against typos silently dropping fields).
- `sample_id` uniqueness enforced across the sheet; duplicates fail parse.
- `reads` paths must be relative — any leading `/` fails the row with an explicit "absolute path not allowed; use `--data-dir` root" error.
- All paths in `reads` are resolved against `--data-dir` and existence-checked before any sample starts. Missing files fail the whole run — no partial execution.
- `kingdom` is normalised to lowercase and matched against the enum. Anything else → parse error naming the offending row.
- `sample_receipt_date`, if present, is parsed as ISO 8601 (`YYYY-MM-DD`); malformed values emit a warning and pass through as the raw string.
- Optional columns are individually optional (empty cell is fine); an omitted column in the header simply means every sample has an implicit empty value for it.
- `--data-dir` must exist and be readable; failure to resolve is a fatal pre-flight error.

**Fan-out:** Nextflow `splitCsv(header: true)` on the samplesheet produces
one `tuple(meta, reads)` per row; `meta` carries `sample_id`, `kingdom`, and
any submitter-supplied optional fields (`sample_info`,
`sample_type`, `sample_receipt_date`, `storage_location`). All stages take
`meta` as a first-class channel key so per-sample context is threaded through
to `COLLATE`, where it is written into `metadata.json` and surfaced in
`report.html`.

**Per-sample output layout:**

```
outdir/
├── INT-2026-0007/
│   ├── organelle_assembly.fasta
│   ├── organelle_annotation.gff
│   ├── barcodes.fasta
│   ├── metadata.json
│   └── report.html
├── INT-2026-0008/
│   └── ...
├── run_manifest.json     # samplesheet snapshot + reference bundle version + pipeline commit
└── run-report.html  # cross-sample summary: N samples, success/no-recovery breakdown
```
