# Task 5 — P1 stages 1 & 4: `NANOPLOT_RAW` and `NANOPLOT_CLEAN`

**Phase:** P1 (first real workflow tool stages — from
[spec §6](../spec/06-phases.md)).
**Goal:** Replace the P0 stubs for stages 1 and 4 with real NanoPlot
invocations that produce read-length / quality metrics before and after
the filtering chain (`CHOPPER` → `FILTLONG`). Both stages run the same
tool with the same input shape and the same output shape, so they are
implemented together.

**Exit criteria:**

- `nextflow run main.nf -profile test` produces, for every sample:
  - `results/<sample_id>/qc/nanoplot_raw/NanoPlot-report.html` +
    `NanoStats.txt` + PNG plots
  - `results/<sample_id>/qc/nanoplot_clean/NanoPlot-report.html` +
    `NanoStats.txt` + PNG plots
  when `publish_intermediates=true`.
- `NanoStats.txt` in each is a real NanoPlot summary (non-empty, contains
  `Mean read length`, `Median read length`, `Mean read quality`, `Number
  of reads`), not a `touch`ed placeholder.
- Both processes run in the pinned biocontainer
  `quay.io/biocontainers/nanoplot:1.47.1--pyhdfd78af_0`
  ([conf/containers.config](../conf/containers.config)) — no host Python,
  no conda.
- CI green — `-profile test` end-to-end + container-coverage check.

**Not in scope:**

- Feeding NanoStats numbers into `report.html` — the report renderer is
  built in P4 ([spec §6a](../spec/06a-reports.md)). This task's contract
  is: the files exist, are real, and are published at the right path.
- Any change to `CHOPPER` / `FILTLONG` / `RECRUIT` or the channel wiring
  between them. `NANOPLOT_CLEAN` continues to consume the FILTLONG
  output tuple with the same shape as P0.
- Trimming / adapter removal — done upstream by Dorado at basecalling
  ([spec §2 stage 1](../spec/02-stages.md#2-stage-detail)).

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- One container per process; container pinned by SHA/tag, never `latest`.
- Every file crossing a process boundary is a Nextflow `Path`.
- No `${params.<file>}` bare interpolation in `script:` blocks.

---

## 1. Architecture

Both `NANOPLOT_RAW` and `NANOPLOT_CLEAN` are thin wrappers around
`NanoPlot`. They differ only in:

| | `NANOPLOT_RAW` | `NANOPLOT_CLEAN` |
|---|---|---|
| Upstream | `PARSE_SAMPLESHEET.out.reads` | `FILTLONG.out.reads` |
| Publish path | `qc/nanoplot_raw/` | `qc/nanoplot_clean/` |
| Output emit | `reports` | `reports` |

Same input shape (`tuple val(meta), path(reads)`), same tool invocation
shape, same output shape. **Do not factor them into a single shared
subworkflow** — the two stages carry different semantic meaning
(baseline vs. post-clean) and downstream consumers (the report, in P4)
reference them by name. Keeping them as two module files matches the
one-file-per-stage convention already established under
[`modules/local/`](../modules/local/).

## 2. NanoPlot command shape

```bash
NanoPlot \
    --threads ${task.cpus} \
    --fastq ${reads} \
    --outdir nanoplot_${variant} \
    --tsv_stats \
    --no_static \
    --format png \
    --title "${meta.sample_id} — ${variant}"
```

Options rationale:

- `--fastq` — accepts `.fastq` or `.fastq.gz` directly; no need to
  decompress.
- `--tsv_stats` — write `NanoStats.txt` in TSV form so downstream
  Python (P4 report) can parse without regex-ing free text.
- `--no_static` — skip static PNG-only pages; the interactive HTML is
  what humans use, PNGs are still produced for the report inclusion via
  `--format png`.
- `--format png` — inline PNG plots (length histogram, quality per
  read, cumulative yield) that the P4 report will embed.
- `--title` — puts the sample id + variant into the HTML title so
  side-by-side comparison of `raw` vs `clean` is unambiguous.
- **No** `--barcoded` — Dorado has already demultiplexed
  ([spec §2 stage 1](../spec/02-stages.md#2-stage-detail)).
- **No** `--minlength` / `--minqual` — filtering happens in CHOPPER /
  FILTLONG; NanoPlot's job is to describe the reads as they are at each
  checkpoint, not to filter them.

Where `${variant}` is `raw` or `clean` per module.

## 3. `modules/local/nanoplot_raw.nf` (replace stub)

```groovy
// Stage 1 — baseline read quality metrics before any filtering.
// Tool: NanoPlot. See spec §2 stage 1.

process NANOPLOT_RAW {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/qc/nanoplot_raw",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("nanoplot_raw/"), emit: reports

    script:
    """
    NanoPlot \\
        --threads ${task.cpus} \\
        --fastq ${reads} \\
        --outdir nanoplot_raw \\
        --tsv_stats \\
        --no_static \\
        --format png \\
        --title "${meta.sample_id} — raw"
    """

    stub:
    """
    mkdir -p nanoplot_raw
    touch nanoplot_raw/NanoPlot-report.html
    touch nanoplot_raw/NanoStats.txt
    """
}
```

## 4. `modules/local/nanoplot_clean.nf` (replace stub)

Identical shape, `clean` substitutions:

```groovy
process NANOPLOT_CLEAN {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/qc/nanoplot_clean",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("nanoplot_clean/"), emit: reports

    script:
    """
    NanoPlot \\
        --threads ${task.cpus} \\
        --fastq ${reads} \\
        --outdir nanoplot_clean \\
        --tsv_stats \\
        --no_static \\
        --format png \\
        --title "${meta.sample_id} — clean"
    """

    stub:
    """
    mkdir -p nanoplot_clean
    touch nanoplot_clean/NanoPlot-report.html
    touch nanoplot_clean/NanoStats.txt
    """
}
```

## 5. `conf/containers.config`

Replace the combined `python:3.12-slim` placeholder:

```groovy
withName: 'NANOPLOT_RAW|NANOPLOT_CLEAN' {
    container = 'quay.io/biocontainers/nanoplot:1.47.1--pyhdfd78af_0'
}
```

Pinned tag, no `latest` ([spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)).

## 6. Resource labels

`process_low` is correct — NanoPlot on a single ONT run (a few Gb of
FASTQ) uses < 2 GB RAM and finishes in 1–3 minutes with 2–4 threads.
No change to [`conf/base.config`](../conf/base.config).

## 7. Tests

### 7.1 Integration — `-profile test`

Add assertions to the existing `-profile test` shell test step (or a
new one under `tests/`):

- Both `nanoplot_raw/NanoStats.txt` and `nanoplot_clean/NanoStats.txt`
  exist under `results/<sample_id>/qc/`.
- Both are non-empty and contain the string `Number of reads` (proves
  NanoPlot actually ran; a `touch`ed stub would fail this).
- Both `NanoPlot-report.html` files exist and are > 10 kB (a real
  NanoPlot HTML is ~1 MB with embedded assets; a stub `touch` gives 0
  bytes).

No unit tests — these processes are pure tool wrappers with no custom
logic ([spec §2.2](../spec/02-stages.md#22-custom-logic-components)).

### 7.2 Manual smoke check

On a real fixture from `tests/data/`:

```bash
nextflow run main.nf -profile test -resume
xdg-open results/TEST-PLANT-01/qc/nanoplot_raw/NanoPlot-report.html
xdg-open results/TEST-PLANT-01/qc/nanoplot_clean/NanoPlot-report.html
```

Sanity-eyeball: read-length distribution and mean Q shift as expected
between raw and clean (fewer reads, higher mean Q, tighter length
distribution after CHOPPER + FILTLONG).

## 8. Deliverables checklist

- [ ] `modules/local/nanoplot_raw.nf` — real `script:` block (stub
      block retained for `-stub-run` mode).
- [ ] `modules/local/nanoplot_clean.nf` — same.
- [ ] `conf/containers.config` — `NANOPLOT_RAW|NANOPLOT_CLEAN` pointed
      at `quay.io/biocontainers/nanoplot:1.47.1--pyhdfd78af_0`.
- [ ] `-profile test` produces real, non-empty NanoStats + HTML for
      both raw and clean, for every sample.
- [ ] CI green — container-coverage check + `-profile test` +
      lint.

## 9. Notes / non-issues

- **`--tsv_stats` file naming.** NanoPlot 1.47.x writes `NanoStats.txt`
  regardless of `--tsv_stats`; the flag only changes the format inside
  the file (TSV vs. free text). Downstream consumers should parse it
  as TSV.
- **Container size.** ~200 MB pulled (NanoPlot + pandas + matplotlib +
  seaborn). Acceptable — the biocontainer is the pinned upstream and
  matches what nf-core modules use for the same tool.
- **`--no_static` interaction with `--format png`.** These are not in
  conflict — `--no_static` disables the static-HTML variant of the
  report page; `--format png` sets the format used for inline plot
  images embedded in the interactive HTML. Both together = interactive
  HTML report with PNG plots (not SVG), which is what we want for
  P4-report embedding.
- **Empty-reads edge case.** If a sample recruits zero reads after
  CHOPPER + FILTLONG (upstream filter too aggressive), NanoPlot exits
  non-zero. This should never happen with the current defaults on real
  ONT data, but if it does surface it is a real signal, not something
  to swallow — leave the default `errorStrategy` so the sample fails
  cleanly rather than silently producing an empty QC bundle. Revisit
  with `errorStrategy 'ignore'` only if P1 test data exposes it.
