# Task 7 — P1 stage 3: `FILTLONG`

**Phase:** P1 (from [spec §6](../spec/06-phases.md)).
**Goal:** Replace the P0 stub for stage 3 with a real `filtlong`
invocation that applies an identity-weighted, percentile-based selection
to the chopper-filtered reads. Purpose: at skim depth, low-identity long
reads swamp the metaFlye assembly graph — filtlong scores each read by
window quality × length and keeps the top N% by that composite score.

**Exit criteria:**

- `nextflow run main.nf -profile test` produces, for every sample:
  `results/<sample_id>/qc/filtlong/<sample_id>.filtlong.fastq.gz` when
  `publish_intermediates=true` — a real gzipped FASTQ, not a `touch`ed
  placeholder.
- Selection parameters (`--min_length`, `--keep_percent`) come from
  [nextflow.config](../nextflow.config#L29-L31) and are visible in the
  process command line via `params` interpolation.
- Process runs in the pinned biocontainer
  `quay.io/biocontainers/filtlong:0.3.1--h077b44d_0`
  ([conf/containers.config](../conf/containers.config#L42-L45)) — no
  host filtlong, no conda.
- The output FASTQ is a valid single-member gzip file that `zcat` reads
  end-to-end without error. Read count is ≤ chopper output count.
- CI green — `-profile test` end-to-end + container-coverage check.

**Not in scope:**

- Reference-guided filtering (`--assembly` / `-a`) — that would require
  an on-hand reference, defeating the species-agnostic design
  ([brief.md §3.1](../brief.md)). Reference-anchored filtering happens
  at `RECRUIT` against the kingdom panel, not here.
- Trimming — done by Dorado upstream ([spec §2 stage 1](../spec/02-stages.md#2-stage-detail)).
- Any change to `CHOPPER` / `NANOPLOT_CLEAN` / `RECRUIT` or the channel
  wiring between them. `NANOPLOT_CLEAN` continues to consume the
  `FILTLONG.out.reads` tuple with the same shape as P0.

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pinned by tag, never `latest`.
- Every file crossing a process boundary is a Nextflow `Path`.
- No `${params.<file>}` bare interpolation in `script:` blocks; scalar
  params (`params.min_read_length`, `params.filtlong_keep_percent`) are
  fine.

---

## 1. `filtlong` command shape

Filtlong reads FASTQ directly (handles `.gz` natively) and writes
FASTQ on stdout:

```bash
filtlong \
    --min_length ${params.min_read_length} \
    --keep_percent ${params.filtlong_keep_percent} \
    ${reads} \
    | gzip > ${meta.sample_id}.filtlong.fastq.gz
```

Options rationale:

- `--min_length N` — hard length floor. Default `500`
  ([nextflow.config §L29](../nextflow.config#L29)); shared with
  chopper's `--minlength` for consistency across the QC chain.
  Duplication is intentional — chopper's cut is on the raw pool,
  filtlong's is on the chopper-cleaned pool; both need the floor
  because filtlong's percentile calc will happily keep the "best" of a
  bad batch if you don't tell it otherwise.
- `--keep_percent 95` — keep the top 95% of reads by filtlong's
  composite score (window mean Q × read length). Default in
  [nextflow.config §L31](../nextflow.config#L31). 95% is a light cut
  that removes clearly-worst reads without depleting yield at skim
  depth. Revisit downwards (e.g. 90%) only if METAFLYE tangles on real
  data ([spec §9](../spec/07-open-questions.md)).
- Positional argument = the input FASTQ path. Filtlong opens it
  directly; no `zcat` needed.
- **No** `--target_bases` — that would clash with `--keep_percent` and
  duplicate the coverage-gate subsampling logic
  ([spec §2.1](../spec/02-stages.md#21-coverage-gate-2-stage-6)).
  Read-count management is `COVERAGE_GATE`'s job.
- **No** `--min_mean_q` — chopper already applied a mean-Q floor
  ([task 6](6_chopper.md)); adding another here is duplicated effort
  and slower (filtlong's per-read scan is more expensive than
  chopper's).
- **No** `-a` / `--assembly` — species-agnostic, see "Not in scope".
- Filtlong is single-threaded — no `--threads` option exists.

## 2. `modules/local/filtlong.nf` (replace stub)

```groovy
// Stage 3 — identity-weighted top-quality read selection.
// Tool: Filtlong. See spec §2 stage 3.

process FILTLONG {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/qc/filtlong",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample_id}.filtlong.fastq.gz"), emit: reads

    script:
    """
    filtlong \\
        --min_length ${params.min_read_length} \\
        --keep_percent ${params.filtlong_keep_percent} \\
        ${reads} \\
        | gzip > ${meta.sample_id}.filtlong.fastq.gz
    """

    stub:
    """
    touch ${meta.sample_id}.filtlong.fastq.gz
    """
}
```

Notes:

- Filtlong keeps the full read set in memory during scoring. For a few
  Gb of ONT FASTQ this is comfortably under the `process_low` 4 GB
  ceiling. On a full nuclear-depth run (rare for this pipeline) it
  might exceed; bump to `process_medium` if that surfaces.
- Do **not** add `set -euo pipefail` — same rationale as
  [task 6 §2](6_chopper.md#2-modulesocalchoppernf-replace-stub) — early
  exit from `gzip` under `pipefail` would false-fail the task.

## 3. `conf/containers.config`

Replace the P0 `python:3.12-slim` placeholder:

```groovy
withName: 'FILTLONG' {
    container = 'quay.io/biocontainers/filtlong:0.3.1--h077b44d_0'
}
```

Pinned tag, no `latest` ([spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)).

## 4. Resource labels

`process_low` is the current label. Filtlong is single-threaded and
uses ~1–2 GB RAM for typical ONT skim runs — well under the 4 GB
`process_low` ceiling in [conf/base.config](../conf/base.config). No
change needed unless a full-depth run surfaces the memory ceiling.

## 5. Tests

### 5.1 Integration — `-profile test`

Add assertions to the existing `-profile test` shell test step:

- `results/<sample_id>/qc/filtlong/<sample_id>.filtlong.fastq.gz`
  exists and is > 0 bytes for every sample.
- `zcat` on the output completes cleanly (exit 0).
- Line count of the output is a positive multiple of 4.
- Line count is ≤ chopper output count and ≤ `keep_percent` fraction
  of it (with slack for the `--min_length` floor). For the synthetic
  test fixtures the numbers are too small to meaningfully assert the
  95% ratio — skip that check under `-profile test` and rely on the
  real-data smoke check for depth accuracy.

No unit tests — filtlong is a pure tool wrapper with no custom logic
([spec §2.2](../spec/02-stages.md#22-custom-logic-components)).

### 5.2 Manual smoke check on real data

After [task 5](5_nanoplot.md) and [task 6](6_chopper.md) land, compare
`NanoStats.txt` in `qc/nanoplot_raw/` vs `qc/nanoplot_clean/`:

- Read count in clean ≈ 95% × chopper output count (slightly less if
  many reads hit the `--min_length` floor).
- Mean quality in clean ≥ mean quality in chopper output (filtlong
  favours high-Q reads).
- Length distribution in clean should be shifted right relative to
  chopper output (filtlong favours long reads).

## 6. Deliverables checklist

- [ ] [`modules/local/filtlong.nf`](../modules/local/filtlong.nf) —
      real `script:` block (stub block retained for `-stub` mode).
- [ ] [`conf/containers.config`](../conf/containers.config) — `FILTLONG`
      pointed at `quay.io/biocontainers/filtlong:0.3.1--h077b44d_0`.
- [ ] `-profile test` produces a real gzipped FASTQ under
      `results/<sample_id>/qc/filtlong/` for every sample when
      `publish_intermediates=true`.
- [ ] CI green — container-coverage check + `-profile test` + lint.

## 7. Notes / non-issues

- **Empty-reads edge case.** Same as [task 6 §7](6_chopper.md#7-notes--non-issues) —
  filtlong on empty input produces empty output cleanly; COVERAGE_GATE
  soft-fails the sample downstream. No special handling here.
- **Why not merge chopper + filtlong into one process.** They are
  distinct enough to keep separate: chopper is hard-threshold, single
  pass, C++; filtlong is percentile, two-pass, Rust. Separate modules
  match the [spec §1 pipeline flow](../spec/01-pipeline-flow.md#1-pipeline-flow),
  make per-stage NanoPlot diagnostics tractable (would need a third
  intermediate stage otherwise), and keep container size / update
  cadence independent.
- **Ordering.** `chopper → filtlong` is deliberate; see
  [task 6 §7](6_chopper.md#7-notes--non-issues) for the rationale.
