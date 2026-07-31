# Task 1 — P0 Scaffold: Nextflow plumbing

**Phase:** P0 (from [plan.md §6](../../plan.md))
**Goal:** Nextflow DSL2 skeleton with every pipeline stage present as a stub
process, channels wired end-to-end, params + samplesheet parsing working,
container plan sketched, CI lint green.

**Exit criteria (from plan.md §6):**
`nextflow run main.nf -profile test` completes a no-op end-to-end pass over
a 2-sample synthetic samplesheet, producing the per-sample output layout
described in [plan.md §0](../../plan.md) with empty/placeholder files.

**Not in scope:** any real tool invocation. All processes emit `touch`ed
placeholder outputs. Real logic arrives in P1+ (see [plan.md §6](../../plan.md)).

**Cross-cutting rules (from [plan.md §1a](../../plan.md)):**

- **Every process runs in a container** — including stubs. No
  host-installed tools, no conda fallbacks. Even a `touch`-only stub
  declares a `container` directive. `-profile test` must pull real
  images and execute stubs inside them; a process without a container
  is a CI failure.
- **File inputs must be stageable.** On remote executors (Azure, AWS
  Batch, K8s, HPC scratch) Nextflow only stages files it sees as a
  `Path` object. Two acceptable patterns
  ([plan.md §1a](../../plan.md)):
  1. **Preferred:** materialise the file into a value channel
     (`ch_x = Channel.value(file(params.x))`) and declare `input: path x`
     on the process. Per-sample data and anything that fans out **must**
     use this pattern.
  2. **Also fine for workflow-wide references** (BLAST DBs, kingdom
     refs, locus panel): inline `${file(params.x)}` in the `script:`
     block — `file()` returns a stageable `Path`, so the executor
     stages it correctly.

  **Bare `${params.x}` in a script block is a bug** — that interpolates
  to a plain string with nothing for the executor to stage. Also no
  `$projectDir`/`$baseDir` file access at runtime. A P0 stub whose
  `input:` block omits a file the real process will need (pattern 1) is
  a wiring bug that must be fixed now, not deferred.

---

## 1. Target repo layout

```
wf5/
├── main.nf                     # entry workflow; wires channels
├── nextflow.config             # profiles, params defaults, container map
├── params.yml                  # example user params (checked in)
├── .env.azure                  # Required env vars for Azure profile
├── conf/
│   ├── base.config             # resource defaults per label
│   ├── test.config             # -profile test (2-sample synthetic run)
│   ├── azure.config            # -profile azure (primary executor for dev/prod)
│   └── containers.config       # one entry per tool image
├── docs/
│   ├── azure/
│   │   └── README.md
│   └── reference-data.md       # Guide to building reference data
├── modules/
│   └── local/
│       ├── parse_samplesheet.nf
│       ├── nanoplot_raw.nf
│       ├── chopper.nf
│       ├── filtlong.nf
│       ├── nanoplot_clean.nf
│       ├── recruit.nf
│       ├── coverage_gate.nf
│       ├── metaflye.nf
│       ├── medaka.nf
│       ├── bandage_ng.nf
│       ├── bin_target.nf
│       ├── blast_validate.nf
│       ├── annotate.nf
│       ├── miniprot_extract.nf
│       ├── ogdraw.nf
│       ├── collate.nf
│       └── run_report.nf
├── scripts/                    # helper scripts invoked from processes
│   ├── tests/
│   │   └── README.md
│   ├── requirements.txt
│   ├── .gitignore              # Python gitignore
│   └── Dockerfile              # Stub for now; `FROM alpine:3.24`
├── assets/
│   └── samplesheet.schema.json # nf-schema JSON schema for samples.csv
├── tests/
│   └── data/
│       ├── samples.csv         # 2 rows: 1 plant, 1 animal, both point at
│       │                       # tiny synthetic fastqs below
│       ├── plant.fastq.gz      # a handful of reads — no real content needed
│       └── animal.fastq.gz
├── wf-report-boilerplate/      # existing; unchanged in this task
├── brief.md                    # existing
├── plan.md                     # existing
├── tests/
│   └── README.md
└── tasks/
    └── 1_scaffold.md           # this file
```

## 2. `main.nf` — workflow shape

Use DSL2. One `workflow {}` block that mirrors the flow diagram in
[plan.md §1](../../plan.md). Key requirements:

1. **Samplesheet fan-out.** `ch_samplesheet | splitCsv(header: true)`
   (with `ch_samplesheet = Channel.value(file(params.samplesheet))`),
   validated against `assets/samplesheet.schema.json` via **nf-schema**
   (`nextflow.enable.dsl=2`, `plugins { id 'nf-schema' }`).
   Emit `tuple(meta, reads_list)` per row where `meta = [sample_id: ...,
   kingdom: ..., sample_info: ..., sample_type: ...,
   sample_receipt_date: ..., storage_location: ...]`.
   Resolve each pipe-delimited `reads` entry against `params.data_dir`;
   reject absolute paths at parse time per [plan.md §0](../../plan.md).
2. **Reads concat.** If a row has >1 fastq, concat them into a single file
   before `NANOPLOT_RAW`. For P0 this can be a `cat` sub-process; keep it
   inside `parse_samplesheet.nf` or a small local process — decide by
   readability.
3. **Linear chain** through `NANOPLOT_RAW` → `CHOPPER` → `FILTLONG` →
   `NANOPLOT_CLEAN` → `RECRUIT` → `COVERAGE_GATE`.
4. **Branch on coverage gate.** `COVERAGE_GATE` emits `tuple(meta, gated_fastq, status_json)`.
   Split with:
   ```groovy
   gated.branch {
       ok:     it[2].text.contains('"status": "ok"')  // for the stub;
                                                       // replace with a
                                                       // proper JsonSlurper
                                                       // read in P1
       failed: true
   }
   ```
   Downstream stages (`METAFLYE` … `OGDRAW`) consume the `ok` branch only.
   Both branches feed into `COLLATE`.
5. **Sequential downstream chain** through `METAFLYE` → `MEDAKA` (opt-in
   via `params.polish`; use `if (params.polish)` around the process call,
   passthrough otherwise) → `BANDAGE_NG` → `BIN_TARGET` → `BLAST_VALIDATE`
   → `ANNOTATE` → `MINIPROT_EXTRACT` → `OGDRAW`.
6. **`COLLATE`** joins the ok- and failed-branch results by `meta.sample_id`
   and emits one per-sample bundle.
7. **`RUN_REPORT`** is the join point: `.collect()` on all `COLLATE`
   outputs, emit run-level `run-report.html` + `run_manifest.json`.

Wire `meta` through every channel as the first tuple element — never lose
it. Cross-sample failure isolation ([plan.md §2.1.4](../../plan.md)) depends on
this.

## 3. Process stub shape

Every module under `modules/local/` follows this template. Real logic
lives in later phases; here we only prove the channel wiring and file
staging work.

```groovy
process CHOPPER {
    tag        "${meta.sample_id}"
    label      'process_low'
    publishDir "${params.outdir}/${meta.sample_id}/qc", mode: 'copy',
               enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample_id}.chopper.fastq.gz"),
        emit: reads

    script:
    """
    # STUB — real implementation in P1
    touch ${meta.sample_id}.chopper.fastq.gz
    """
}
```

Rules for stubs:

- Always emit **exactly the file shapes** the real process will emit
  (name, extension, tuple layout) so downstream wiring is correct now.
- `touch` files rather than generating fake content — no risk of a later
  reader mistaking placeholder bytes for real data.
- Use `tag "${meta.sample_id}"` on every per-sample process so
  `nextflow log` is readable.
- Use `label` (e.g. `process_low` / `process_medium` / `process_high`)
  with resources defined once in `conf/base.config`.
- `container` directive resolved centrally by `conf/containers.config`
  keyed off process name — **mandatory for every process** ([plan.md §1a](../../plan.md)),
  stubs included. Pin every tag; no `latest`. For processes whose real
  tool is not yet wired (most of them in P0), use a small, always-safe
  base image (e.g. `quay.io/biocontainers/coreutils:9.5`) so `touch`
  works; replace with the tool image as each stage lands in P1+.
- `errorStrategy 'ignore'` set **only** on `COVERAGE_GATE`
  (per [plan.md §2.1.4](../../plan.md)) — nowhere else.

### Per-process I/O shapes (all stubs)

| # | Process | Input tuple | Output tuple(s) |
|---|---|---|---|
| 0 | `PARSE_SAMPLESHEET` | (channel from samplesheet) | `tuple(meta, path(reads))` — one item per row, reads concatenated if multi-file |
| 1 | `NANOPLOT_RAW` | `tuple(meta, reads)` | `tuple(meta, path("nanoplot_raw/"))` |
| 2 | `CHOPPER` | `tuple(meta, reads)` | `tuple(meta, path("*.chopper.fastq.gz"))` |
| 3 | `FILTLONG` | `tuple(meta, reads)` | `tuple(meta, path("*.filtlong.fastq.gz"))` |
| 4 | `NANOPLOT_CLEAN` | `tuple(meta, reads)` | `tuple(meta, path("nanoplot_clean/"))` |
| 5 | `RECRUIT` | `tuple(meta, reads)` + `path(kingdom_refs)` | `tuple(meta, path("*.recruited.fastq.gz"))` |
| 6 | `COVERAGE_GATE` | `tuple(meta, reads)` | `tuple(meta, path("*.gated.fastq.gz"), path("sample_status.json"), path("coverage.json"))` |
| 7 | `METAFLYE` | `tuple(meta, reads)` | `tuple(meta, path("assembly.fasta"), path("assembly_graph.gfa"), path("assembly_info.txt"))` |
| 8 | `MEDAKA` | `tuple(meta, assembly, reads)` | `tuple(meta, path("*.polished.fasta"))` |
| 9 | `BANDAGE_NG` | `tuple(meta, gfa)` | `tuple(meta, path("*.graph.png"))` |
| 10 | `BIN_TARGET` | `tuple(meta, assembly, gfa, info)` + `path(kingdom_refs)` | `tuple(meta, path("target.fasta"), path("secondaries.tsv"), optional path("path[12].fasta") for plant)` |
| 11 | `BLAST_VALIDATE` | `tuple(meta, target_fasta)` + `path(blast_db)` | `tuple(meta, path("*.blast.tsv"))` |
| 12 | `ANNOTATE` | `tuple(meta, target_fasta)` | `tuple(meta, path("*.gff"), path("*.gbk"))` |
| 13 | `MINIPROT_EXTRACT` | `tuple(meta, target_fasta)` + `path(locus_panel)` | `tuple(meta, path("barcodes.fasta"), path("*.coords.gff"), path("*.validation.tsv"))` |
| 14 | `OGDRAW` | `tuple(meta, gbk)` | `tuple(meta, path("*.map.svg"))` |
| 15 | `COLLATE` | *joined per-sample bundle* | `tuple(meta, path("${meta.sample_id}/"))` — the whole outdir layout from [plan.md §0](../../plan.md) |
| 16 | `RUN_REPORT` | `.collect()` of all `COLLATE` outputs + `path(samplesheet)` | `path("run-report.html"), path("run_manifest.json")` |

**Value-channel plumbing:** at the top of `main.nf`, materialise each
reference param into a value channel exactly once, then reuse:

```groovy
ch_kingdom_refs = Channel.value(file(params.kingdom_refs))
ch_locus_panel  = Channel.value(file(params.locus_panel))
ch_blast_db     = Channel.value(file(params.blast_db))
ch_samplesheet  = Channel.value(file(params.samplesheet))
```

Passing `file(params.x)` guarantees the file is a stageable Nextflow
`Path` object — the executor knows to copy it to the compute node.
Never reference `params.kingdom_refs` (or any other file param)
directly inside a process `script:` block.

For `COVERAGE_GATE` the stub script writes a hardcoded
`{"status": "ok"}` JSON so the ok-branch flows end-to-end. Leave a
`TODO` comment noting P1 will make this a real decision, and consider
adding a second `-profile test_gated_fail` in `conf/` that writes
`{"status": "low_coverage"}` to prove the failed branch reaches
`COLLATE`.

## 4. Params schema

`nextflow.config` `params` block, all with sensible defaults:

| Param | Type | Default | Notes |
|---|---|---|---|
| `samplesheet` | path | — | required |
| `data_dir` | path | — | required; resolves relative `reads` paths ([plan.md §0](../../plan.md)) |
| `outdir` | path | `./results` | |
| `polish` | bool | `false` | opt-in Medaka ([plan.md §2 stage 8](../../plan.md)) |
| `publish_intermediates` | bool | `false` | if true, intermediate stage outputs are copied under `outdir/<sample_id>/<stage>/` |
| `kingdom_refs` | path | *placeholder* | reference bundle root ([plan.md §4](../../plan.md)); P0 accepts any dir |
| `locus_panel` | path | *placeholder* | Taxodactyl accepted-loci config ([plan.md §4.3](../../plan.md)) |
| `blast_db` | path | *placeholder* | RefSeq organelle DB ([plan.md §4.2](../../plan.md)) |
| `min_read_length` | int | — | CHOPPER threshold (P1 wires) |
| `min_mean_q` | int | — | CHOPPER threshold (P1 wires) |
| `filtlong_keep_percent` | int | 95 | FILTLONG param ([plan.md §2 stage 3](../../plan.md)) |

Ship `assets/samplesheet.schema.json` matching the samplesheet contract
in [plan.md §0](../../plan.md): `sample_id` pattern `^[A-Za-z0-9_.-]+$`,
`kingdom` enum `[plant, animal]`, `reads` required, optional columns
allowed. nf-schema will validate at workflow entry and produce readable
errors.

## 5. `conf/` details

- **`base.config`** — resource labels only:
  ```groovy
  process {
      withLabel: process_low    { cpus = 2;  memory = 4.GB }
      withLabel: process_medium { cpus = 4;  memory = 16.GB }
      withLabel: process_high   { cpus = 8;  memory = 32.GB }
  }
  ```
- **`containers.config`** — one entry per process, keyed by process name.
  **Mandatory for every process in P0** ([plan.md §1a](../../plan.md)) —
  stubs point at a small placeholder image (e.g.
  `quay.io/biocontainers/coreutils:9.5`); real tool images replace them
  as each stage is wired in P1+. Pin every tag; no `latest`. A missing
  entry is a CI failure (see §7). **Do not bake custom Python source
  into the placeholder images**: the "deps in image, code at runtime"
  pattern ([plan.md §2.2](../../plan.md)) applies from P0 onward — helper
  scripts under `bin/` are auto-staged by Nextflow onto the container
  `PATH` at execution time. `images.yml` should therefore rebuild only
  when `scripts/requirements.txt` or a Dockerfile changes
  ([plan.md §5b](../../plan.md)), not on every Python edit.
- **`test.config`** — points `samplesheet`, `data_dir`, `kingdom_refs`,
  `locus_panel`, `blast_db` at the fixtures under `tests/data/`.
  Enables the container runtime (Docker by default; Singularity/Apptainer
  optional via a sibling `test_singularity` profile). **Do not disable
  containers** — even stubs execute inside their declared image.

## 6. Test fixtures

- `tests/data/samples.csv` — 2 rows exercising both kingdoms and the
  optional columns:
  ```
  sample_id,kingdom,reads,sample_info,sample_type,sample_receipt_date,storage_location
  TEST-PLANT-01,plant,plant.fastq.gz,synthetic plant test,leaf,2026-07-24,
  TEST-ANIMAL-01,animal,animal.fastq.gz,synthetic animal test,whole,2026-07-24,
  ```
- `tests/data/plant.fastq.gz`, `animal.fastq.gz` — a handful of dummy
  reads each (single-record `@stub` entry is fine; content is never
  consumed by stubs).

## 7. CI

GitHub Actions workflow `.github/workflows/lint.yml`:

- `nf-core lint` (or bare `nextflow lint main.nf` if we skip the
  nf-core scaffolding) — catches syntax + convention issues.
- **Container-coverage check.** A small script (or lint step) asserts
  every process in `modules/local/` resolves to a `container` directive
  via `containers.config`, and that no tag is `latest`
  ([plan.md §1a](../../plan.md)). A process without a container fails CI.
- `nextflow run main.nf -profile test -stub` — runs the stub end-to-end
  **with the container runtime enabled**, so image-pull failures surface
  in CI. All processes have `stub:` blocks identical to their `script:`
  block in P0, so `-stub` is trivially green. Once real logic lands
  per-stage, the `script:` will grow but `stub:` stays as a placeholder
  for CI.

## 8. Deliverables checklist

- [ ] Repo layout above created; no directories empty (use `.gitkeep`).
- [ ] `main.nf` wires all 17 processes with `meta` threaded through.
- [ ] Coverage-gate branch split working end-to-end (both `ok` and
      `failed` reach `COLLATE`).
- [ ] `-profile test` completes successfully; `results/TEST-PLANT-01/`
      and `results/TEST-ANIMAL-01/` contain the per-sample layout from
      [plan.md §0](../../plan.md) (populated with placeholder files);
      `results/run-report.html` and `results/run_manifest.json` exist.
- [ ] nf-schema validation rejects an invalid samplesheet (e.g. bad
      kingdom, absolute path in `reads`) with a clear error.
- [ ] Every process in `modules/local/` resolves to a pinned container
      image via `containers.config` ([plan.md §1a](../../plan.md)); CI
      container-coverage check passes.
- [ ] `-profile test` runs stubs **inside their containers** (Docker
      enabled), not on the host.
- [ ] No process `script:` block contains bare `${params.<file>}`
      interpolation ([plan.md §1a](../../plan.md) staging rule). CI
      lint: `grep -rnE '\$\{?params\.' modules/local/` and fail on any
      match that is not wrapped in `file(...)` or gated behind an
      allow-list of scalar params (`params.polish`, thread counts,
      etc.). `input: path` declarations and `${file(params.x)}` inline
      references are both fine.
- [ ] CI green on push.

## 9. Open questions to resolve before starting

1. **nf-core scaffold, or bare Nextflow?** nf-core buys us a lot of
   scaffolding (lint, schema, module conventions, test harness) at the
   cost of some ceremony. Recommend **nf-core** given the number of
   modules and the CI story, but flag if the team prefers to keep it
   lean.
2. **Container runtime.** Docker vs Singularity/Apptainer for local dev
   and for the HPC target. Containers are **mandatory from P0**
   ([plan.md §1a](../../plan.md)), so this must be resolved before
   scaffolding starts. Recommend Docker as the default profile
   (matches CI) with a `test_singularity` sibling profile for HPC-adjacent
   dev machines.
3. **`MEDAKA` passthrough.** When `params.polish = false` (default),
   should the workflow use an `if` around the process call
   (assembly channel flows straight to `BANDAGE_NG`), or always call
   `MEDAKA` with a no-op stub that just re-emits input? Recommend the
   `if` — cheaper, clearer in `-log`.
