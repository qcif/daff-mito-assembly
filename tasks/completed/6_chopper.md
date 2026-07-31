# Task 6 — P1 stage 2: `CHOPPER`

**Phase:** P1 (from [spec §6](../../spec/06-phases.md)).
**Goal:** Replace the P0 stub for stage 2 with a real `chopper` invocation
that applies a minimal length + mean-quality filter to the raw reads.
Purpose: strip clearly-unusable reads before FILTLONG's percentile-based
selection so FILTLONG's "keep top 95%" operates on a clean pool.

**Exit criteria:**

- `nextflow run main.nf -profile test` produces, for every sample:
  `results/<sample_id>/qc/chopper/<sample_id>.chopper.fastq.gz` when
  `publish_intermediates=true` — a real gzipped FASTQ, not a `touch`ed
  placeholder.
- Filter parameters (`--min_read_length`, `--min_mean_q`) come from
  [nextflow.config](../../nextflow.config#L29-L30) and are visible in the
  process command line via [`params`](../../nextflow.config) interpolation
  (not hardcoded in the module).
- Process runs in the pinned biocontainer
  `quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0`
  ([conf/containers.config](../../conf/containers.config#L38-L41)) — no
  host chopper, no conda.
- The output FASTQ is a valid single-member gzip file that `zcat` reads
  end-to-end without error. Read count is ≤ input read count.
- CI green — `-profile test` end-to-end + container-coverage check.

**Not in scope:**

- Adapter / barcode trimming — done upstream by Dorado at basecalling
  time ([spec §2 stage 1](../../spec/02-stages.md#2-stage-detail)). Do
  **not** pass `--headcrop` / `--tailcrop`.
- Contamination filtering — positive recruitment happens at `RECRUIT`
  ([spec §2 stage 5](../../spec/02-stages.md#2-stage-detail)). Do not pass
  `--contam`.
- Maximum-length filter (`--maxlength`) — organelle reads span the full
  ONT length distribution; capping length would discard genuine target
  reads.
- Any change to `NANOPLOT_RAW` / `FILTLONG` / `NANOPLOT_CLEAN` or the
  channel wiring between them. `FILTLONG` continues to consume the
  `CHOPPER.out.reads` tuple with the same shape as P0.

**Cross-cutting rules (from [spec §1a](../../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pinned by tag, never `latest`.
- Every file crossing a process boundary is a Nextflow `Path`.
- No `${params.<file>}` bare interpolation in `script:` blocks; scalar
  params (`params.min_read_length`, `params.min_mean_q`) are fine.

---

## 1. `chopper` command shape

`chopper` reads FASTQ on stdin and writes filtered FASTQ on stdout:

```bash
zcat ${reads} \
    | chopper \
        --threads ${task.cpus} \
        --quality ${params.min_mean_q} \
        --minlength ${params.min_read_length} \
    | gzip > ${meta.sample_id}.chopper.fastq.gz
```

Options rationale:

- `--quality N` — drops reads with mean Q < N. Default in
  [nextflow.config](../../nextflow.config#L30) is `10`; Dorado SUP output
  is typically ~Q20+ so this is a floor, not a filter.
- `--minlength N` — drops reads shorter than N. Default `500`
  ([nextflow.config §L29](../../nextflow.config#L29)); organelle reads
  worth recruiting are almost always longer than this.
- `--threads` — chopper is IO-bound at low thread counts but scales to
  4+ cheaply.
- **No** `--maxlength` — see "Not in scope".
- **No** `--headcrop` / `--tailcrop` — Dorado already trimmed adapters
  ([spec §2 stage 1](../../spec/02-stages.md#2-stage-detail)).
- **No** `--contam` — positive recruitment does contamination handling
  downstream ([spec §2 stage 5](../../spec/02-stages.md#2-stage-detail)).

The `zcat | chopper | gzip` pattern is idiomatic for chopper — the tool
does not have a `--input`/`--output` file API; it is a stdin/stdout
filter by design.

## 2. `modules/local/chopper.nf` (replace stub)

```groovy
// Stage 2 — length/quality filter.
// Tool: chopper. See spec §2 stage 2.

process CHOPPER {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/qc/chopper",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample_id}.chopper.fastq.gz"), emit: reads

    script:
    """
    zcat ${reads} \\
        | chopper \\
            --threads ${task.cpus} \\
            --quality ${params.min_mean_q} \\
            --minlength ${params.min_read_length} \\
        | gzip > ${meta.sample_id}.chopper.fastq.gz
    """

    stub:
    """
    touch ${meta.sample_id}.chopper.fastq.gz
    """
}
```

Notes:

- `zcat` handles both `.fastq` and `.fastq.gz` inputs — for uncompressed
  input `zcat` just passes bytes through. If the biocontainer omits
  `zcat`, swap to `gunzip -c` (bash builtin — available in every
  biocontainer via BusyBox / coreutils base).
- Do **not** add `set -euo pipefail` to the script block — Nextflow
  already wraps script blocks with `-e` (`errexit`) semantics, and a
  broken pipe from `gzip` on early exit would false-fail the task under
  `pipefail`. If a real broken-pipe scenario surfaces, revisit.

## 3. `conf/containers.config`

Replace the P0 `python:3.12-slim` placeholder:

```groovy
withName: 'CHOPPER' {
    container = 'quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0'
}
```

Pinned tag, no `latest` ([spec §1a](../../spec/01-pipeline-flow.md#1a-engineering-constraints)).

## 4. Resource labels

`process_low` is correct — chopper on a full ONT run (a few Gb of
FASTQ) uses < 1 GB RAM and finishes in 2–5 minutes with 2–4 threads.
No change to [`conf/base.config`](../../conf/base.config).

## 5. Tests

### 5.1 Integration — `-profile test`

Add assertions to the existing `-profile test` shell test step:

- `results/<sample_id>/qc/chopper/<sample_id>.chopper.fastq.gz` exists
  and is > 0 bytes for every sample.
- `zcat` on the output completes cleanly (exit 0) — proves valid gzip.
- Line count of the output is a positive multiple of 4 (FASTQ record =
  4 lines).
- Line count of the output is ≤ line count of the input — chopper
  never adds reads, only drops.

No unit tests — chopper is a pure tool wrapper with no custom logic
([spec §2.2](../../spec/02-stages.md#22-custom-logic-components)).

### 5.2 Manual smoke check on real data

Compare `NanoStats.txt` before and after (once
[task 5](5_nanoplot.md) lands):

- Mean read length should be ≥ raw mean (short reads dropped).
- Mean quality should be ≥ raw mean (low-Q reads dropped).
- Read count should be < raw count.

## 6. Deliverables checklist

- [ ] [`modules/local/chopper.nf`](../../modules/local/chopper.nf) — real
      `script:` block (stub block retained for `-stub` mode).
- [ ] [`conf/containers.config`](../../conf/containers.config) — `CHOPPER`
      pointed at `quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0`.
- [ ] `-profile test` produces a real gzipped FASTQ under
      `results/<sample_id>/qc/chopper/` for every sample when
      `publish_intermediates=true`.
- [ ] CI green — container-coverage check + `-profile test` + lint.

## 7. Notes / non-issues

- **Empty-reads edge case.** If chopper drops all reads (input was
  entirely below thresholds), it exits 0 with an empty output. FILTLONG
  downstream will do the same, then COVERAGE_GATE will soft-fail the
  sample. That is the correct behaviour — a genuinely unusable sample
  becomes a low-coverage soft-fail with clear reporting
  ([spec §2.1](../../spec/02-stages.md#21-coverage-gate-2-stage-6)),
  not a pipeline error.
- **Why chopper before filtlong, not the reverse.** Chopper applies a
  hard threshold (drop below Q10 / 500 bp). Filtlong keeps the top 95%
  of what remains, weighted by identity. Reversing the order would let
  filtlong's percentile calculation be dragged down by clearly-junk
  reads that chopper would have removed for free. Order matches
  [spec §1 pipeline flow](../../spec/01-pipeline-flow.md#1-pipeline-flow).
- **Tool spelling.** The binary is `chopper` (lowercase). No `Chopper`,
  no `nanofilt` — that is a different tool and has been superseded by
  chopper for ONT filtering.
