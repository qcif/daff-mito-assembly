# Task 8 — P1 stage 5: `RECRUIT`

**Phase:** P1 (from [spec §6](../spec/06-phases.md)).
**Goal:** Replace the P0 stub for stage 5 with a real positive-recruitment
step: align QC-passed reads against the sample's target organelle
reference with minimap2, keep the mapped read IDs, and pull those reads
back out of the input FASTQ with seqtk. Reads that don't recruit are
discarded ([brief.md §3.3](../brief.md)).

Per [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)
each sample-sheet row targets exactly one organelle
(`meta.assembly_target ∈ {animal_mt, plant_pt, plant_mt}`), so RECRUIT
aligns against a single index — no per-organelle fan-out inside a
sample.

**Prerequisite:** the `assembly_target` migration ([task 9](9_assembly_target_migration.md))
must land first so `meta.assembly_target` is populated on the channel.

**Exit criteria:**

- `nextflow run main.nf -profile test` produces, for every sample:
  `results/<sample_id>/recruit/<sample_id>.recruited.fastq.gz` when
  `publish_intermediates=true` — a real gzipped FASTQ, not a `touch`ed
  placeholder.
- Recruited-read counts for the CI fixtures are non-zero (test refs
  built from GetOrganelleDB subsets that include the source organism —
  see §5).
- Recruitment aligns against exactly one index per run, selected by
  `meta.assembly_target`: `animal_mt` → `animal_mt.mmi`, `plant_pt` →
  `plant_pt.mmi`, `plant_mt` → `plant_mt.mmi`.
- Process runs in the mulled biocontainer
  `quay.io/biocontainers/mulled-v2-9278ceb357570ba6e25522f5a16e2a4d3ba61a68:74b9019c6ae38f81f95dd09e981151bb2d1028ee-0`
  ([conf/containers.config](../conf/containers.config)) — already
  pinned. No host tools.
- Output FASTQ preserves original read names, sequences, and quality
  strings (proof: `seqtk subseq` copies bytes from the input; it does
  not round-trip through samtools reconstruction). Read count ≤ input
  filtlong count.
- CI green — `-profile test` end-to-end + container-coverage check.

**Not in scope:**

- The stricter PAF-based ptGAUL filter ([spec §2 stage 5](../spec/02-stages.md#2-stage-detail)
  "Optional stricter filter") — deferred to a P1 benchmark task.
- Coverage estimation / subsampling / soft-fail — that is COVERAGE_GATE's
  job ([spec §2.1](../spec/02-stages.md#21-coverage-gate-2-stage-6)),
  landing in a separate task.
- Reference bundle building — that is [task 3](3_refdata.md). This task
  consumes an already-built bundle at `params.kingdom_refs`.
- Any change to BIN_TARGET, which also consumes `ch_kingdom_refs` — its
  own real-implementation task will handle its per-kingdom logic.

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pinned by SHA, never `latest`.
- Reference bundle enters the process via an `input: path` channel —
  never bare `${params.kingdom_refs}` in the script block.
- Every file crossing a process boundary is a Nextflow `Path`.

---

## 1. `params.organelle_refs` semantic change

Currently `params.kingdom_refs` (renamed by [task 9](9_assembly_target_migration.md)
to `params.organelle_refs`) points at a single `.mmi` file (P0
placeholder). This task changes it to point at a **directory** matching
the bundle layout defined in [task 3 §1](3_refdata.md#1-target-layout):

```
<params.organelle_refs>/
├── animal_mt.mmi        # used when meta.assembly_target == 'animal_mt'
├── plant_pt.mmi         # used when meta.assembly_target == 'plant_pt'
└── plant_mt.mmi         # used when meta.assembly_target == 'plant_mt'
```

Filename := `<assembly_target>.mmi` — 1:1 mapping.

`Channel.value(file(params.organelle_refs))` in `main.nf` handles
directories transparently — Nextflow stages the directory as a symlink
into the work dir. The downstream `BIN_TARGET(..., ch_organelle_refs)`
call also gets the directory, which it will consume the same way when
its real-implementation task lands.

## 2. Recruitment command shape

Two-step: (1) align + collect mapped read IDs; (2) `seqtk subseq` back
against the original FASTQ. Rationale: `seqtk subseq` preserves the
original read record byte-for-byte, whereas `samtools fastq` reformats
the record from BAM fields and can drop / munge auxiliary fields.

```bash
minimap2 -ax map-ont -t $T $refs/${assembly_target}.mmi $reads \
    | samtools view -F 4 -q 1 -@ $T \
    | cut -f1 \
    | sort -u > ids.txt

seqtk subseq $reads ids.txt | gzip > $out
```

One index per run — no per-organelle fan-out inside a sample. Plant
samples requiring both plastid and mitogenome are two rows in the
sample sheet ([spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)),
and each row hits this script once with its own `assembly_target`.

Flag rationale:

- `minimap2 -ax map-ont` — SAM output, ONT preset (long, error-prone).
  Spec-mandated at [§2 stage 5](../spec/02-stages.md#2-stage-detail).
- `samtools view -F 4` — drop unmapped reads (flag 4 = unmapped).
- `samtools view -q 1` — drop MAPQ 0 (minimap2 emits MAPQ 0 for
  multi-mappers with no primary; keeping them would let stray host /
  contaminant reads through when they happen to match a conserved
  organelle stretch). Spec-mandated.
- `cut -f1 | sort -u` — deduplicate on read name (multi-mapping reads
  emit multiple alignment records).
- `seqtk subseq` — pulls the *original* FASTQ record for each ID, so
  downstream stages (COVERAGE_GATE + METAFLYE) see unmodified reads.

Empty `ids.txt` → empty `seqtk subseq` output → 20-byte empty gzip →
COVERAGE_GATE soft-fails the sample with `status: "low_coverage"`. That
is the correct behaviour — a sample that recruits zero organelle reads
should not proceed to assembly.

## 3. `modules/local/recruit.nf` (replace stub)

```groovy
// Stage 5 — positive recruitment against kingdom organelle reference panel.
// Tools: minimap2 + samtools + seqtk (mulled biocontainer).
// See spec §2 stage 5 and brief §3.2–3.3.

process RECRUIT {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/recruit",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)
    path organelle_refs

    output:
    tuple val(meta), path("${meta.sample_id}.recruited.fastq.gz"), emit: reads

    stub:
    """
    touch ${meta.sample_id}.recruited.fastq.gz
    """

    script:
    // Target-driven index selection. The sample's assembly_target names
    // exactly one .mmi in the refs bundle.
    """
    minimap2 -ax map-ont -t ${task.cpus} \\
             ${organelle_refs}/${meta.assembly_target}.mmi ${reads} \\
        | samtools view -F 4 -q 1 -@ ${task.cpus} \\
        | cut -f1 \\
        | sort -u > ids.txt

    seqtk subseq ${reads} ids.txt \\
        | gzip > ${meta.sample_id}.recruited.fastq.gz
    """
}
```

Notes:

- Do **not** add `set -euo pipefail` — same rationale as
  [task 6 §2](6_chopper.md#2-modulesocalchoppernf-replace-stub) — a
  broken pipe under `pipefail` on the last stage would false-fail the
  task.
- `process_medium` label is correct — minimap2 uses ~1–4 GB RAM on the
  organelle-scale references and the task benefits from 4–8 threads.
- Ensure the [task 9](9_assembly_target_migration.md) migration is in
  before this task lands — `meta.assembly_target` must exist on the
  channel.

## 4. `conf/test.config` and `conf/test_bad_samplesheet.config`

Both currently point at a single file:

```groovy
kingdom_refs = "${projectDir}/tests/data/refs/kingdom_refs.mmi"
```

Change to the recruit directory:

```groovy
kingdom_refs = "${projectDir}/tests/data/refs/recruit"
```

Delete the empty placeholder `tests/data/refs/kingdom_refs.mmi`
alongside this change.

## 5. Test reference fixtures

The CI fixtures ([tests/data/](../tests/data/)) come from real
organisms — pea aphid (SRR8306868) for animal, Datura stramonium
(SRR11315861, once its download completes) for plant. To get non-zero
recruitment on those fixtures we need small kingdom refs that include
the source organism's organelle sequences.

**Test refs layout:**

```
tests/data/refs/recruit/
├── animal_mt.mmi      # includes A. pisum mt (or a near-relative)
├── plant_pt.mmi       # includes D. stramonium cp (or a near Solanaceae)
├── plant_mt.mmi       # includes D. stramonium mt (or a near Solanaceae)
└── README.md          # source list + regen recipe
```

**Build recipe** (documented in `tests/data/refs/recruit/README.md`):

```bash
# Animal: subset GetOrganelleDB animal_mt.fasta for aphid + close relatives.
docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v /path/to/getorganelle-db/0.0.1/SeedDatabase:/db:ro \
    -v "$PWD/tests/data/refs/recruit":/out \
    quay.io/biocontainers/seqkit:2.10.1--h9ee0642_0 \
    seqkit grep -r -n -p 'aphid|hemiptera|acyrthosiphon' \
        /db/animal_mt.fasta -o /out/animal_mt.fa

# Plant plastid: subset embplant_pt.fasta for Solanaceae.
# ... similar for plant_mt using embplant_mt.fasta ...

# Build minimap2 indices (mulled RECRUIT container has minimap2).
docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v "$PWD/tests/data/refs/recruit":/w -w /w \
    quay.io/biocontainers/mulled-v2-9278ceb357570ba6e25522f5a16e2a4d3ba61a68:74b9019c6ae38f81f95dd09e981151bb2d1028ee-0 \
    bash -c 'for k in animal_mt plant_pt plant_mt; do
                minimap2 -d ${k}.mmi ${k}.fa
             done'
```

Size expectation: subsetted FASTAs are ~100 kB–2 MB each; `.mmi`
indices are similar size. Total addition to git: ~2–6 MB. Acceptable.

**Fallback if source-organism entries are missing from GetOrganelleDB:**
fall back to whole-kingdom subsets (top 20 largest sequences from each
seed FASTA). Fixture reads may not recruit, and COVERAGE_GATE will
soft-fail — that is still a legitimate CI signal that RECRUIT itself
runs cleanly on real data.

## 6. Verification

### 6.1 Integration — `-profile test`

Assertions to add to the `-profile test` shell test step:

- `results/<sample_id>/recruit/<sample_id>.recruited.fastq.gz` exists
  and is > 20 bytes (i.e. non-empty gzip) for at least the animal
  samples (guaranteed by A. pisum mt being in the animal test ref).
- `zcat` on each output completes cleanly.
- Recruited read count ≤ filtlong output count for the same sample.
- Recruited read names are a subset of filtlong output read names
  (proves `seqtk subseq` pulled from the correct input, not somehow
  reconstructed).

### 6.2 Manual smoke check

For plant samples once the real plant fixture lands: submit both a
`plant_pt` and a `plant_mt` row sharing the same reads and confirm
each produces its own recruitment output under a distinct `sample_id`,
with different recruited-read counts reflecting the two organelle
references.

## 7. Deliverables checklist

- [x] [`modules/local/recruit.nf`](../modules/local/recruit.nf) — real
      `script:` block; single minimap2 pass against
      `${organelle_refs}/${meta.assembly_target}.mmi` (stub retained).
- [x] [`conf/test.config`](../conf/test.config) — `organelle_refs` points
      at `tests/data/refs/recruit` (directory).
- [x] [`conf/test_bad_samplesheet.config`](../conf/test_bad_samplesheet.config)
      — same update.
- [x] `tests/data/refs/organelle_refs.mmi` — deleted (renamed empty
      placeholder from task 9 migration, superseded by the recruit dir).
- [x] `tests/data/refs/recruit/{animal_mt,plant_pt,plant_mt}.{fa,mmi}` —
      built from GetOrganelleDB subsets (11 Hemiptera seqs, 10 diverse
      chloroplasts, 1 plant mt). See
      [`tests/data/refs/recruit/README.md`](../tests/data/refs/recruit/README.md).
- [x] `tests/data/refs/recruit/README.md` — sources + regen recipe.
- [x] `-profile test` produces non-empty recruited FASTQ for all 4
      samples (see §9 below for counts).
- [x] CI green — `-profile test` 62/62.

## 8. Notes / non-issues

- **Empty output is a valid state.** A sample whose reads don't recruit
  produces an empty gzip; COVERAGE_GATE soft-fails it in the next
  stage. No `errorStrategy` change here — RECRUIT itself exits 0
  when the tools succeed but return no hits.
- **`-q 1` vs no MAPQ filter.** Spec §2 stage 5 pins `-q 1`. In
  practice minimap2 assigns MAPQ 0 to fully-multi-mapped reads with no
  best primary. Dropping those trims noise from conserved-region hits
  in unrelated organisms; keeping them would let more marginal recruits
  through, which we do not want at the recruitment stage (COVERAGE_GATE
  and BIN_TARGET are downstream cleanups but shouldn't have to
  compensate for lax recruitment).
- **Why `seqtk subseq` not `samtools fastq`.** `samtools fastq`
  reconstructs the FASTQ record from BAM fields — quality strings can
  be reversed if the read aligned to the reverse strand, and any
  auxiliary tags are dropped. `seqtk subseq` copies the original record
  verbatim, which is what downstream QC / assembly expects.
- **No within-sample fan-out.** A plant operator wanting both `plant_pt`
  and `plant_mt` bundles submits two sample-sheet rows sharing the same
  reads (per [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)),
  and each recruits against its own single index. Keeps process count,
  work dirs, and per-target parameters (coverage limits, genome-size
  hint, BLAST DB, genetic code) unambiguous.
- **Container `--user` propagation.** Already handled via
  [nextflow.config §60](../nextflow.config#L60) — no per-process
  config needed.

---

## 9. Findings — recruitment yield on CI fixtures

**The original random-subsample fixtures produced zero recruited reads.**
This is expected: with 70 reads randomly sampled from a 892k-read WGS
run, the probability of drawing a read that maps to the organelle is
≈ 0 at the scales we used. Fixtures were regenerated as
**pre-recruited subsets** — recruit against the test refs from the full
staging FASTQ, then subsample from that organelle-enriched pool.

| Sample | Target | Recruited reads | Source pool size |
|---|---|---|---|
| TEST-ANIMAL-01 | `animal_mt` | **56 / 70** (80%) | 1,242 total in pool |
| TEST-ANIMAL-02 | `animal_mt` | **86 / 105** (82%) | 1,242 (animal + animal_b concat) |
| TEST-PLANT-01-pt | `plant_pt` | **25 / 15** (plant_pt reads) | 53,641 total in pool |
| TEST-PLANT-01-mt | `plant_mt` | **18 / 15** (plant_mt reads) | 14,981 total in pool |

Plant recruited counts exceed the 15-read input because the `plant.fastq.gz`
fixture is a mix of 15 plant_pt + 15 plant_mt reads — some plant_mt reads
also hit the plant_pt reference (shared conserved genes), and vice versa.

**Implications for downstream stages:**

- COVERAGE_GATE will soft-fail all CI fixtures — recruited bases are far
  below the 30× minimum for any target (e.g., 56 reads × ~5 kb/read =
  280 kb vs 17 kb animal mt nominal → ~16× coverage, but these are
  quality-filtered reads so actual coverage is less). At ~56 reads this
  is right on the margin; COVERAGE_GATE behaviour in CI depends on
  exactly what coverage.py computes against the nominal organelle size.
  When COVERAGE_GATE is implemented (task 10), confirm whether the CI
  fixture recruits enough to pass the gate or intentionally soft-fails.
- Real production samples (Tier 2 integration) will see much higher
  recruitment yield. Skim WGS at ~2× nuclear depth gives ~100-200×
  organelle coverage, well above the MIN threshold.
- If CI tests want to exercise the METAFLYE → BIN_TARGET → EXTRACT path
  (not just the soft-fail path), the fixture will need to be extended or
  the test profile will need to override the MIN coverage threshold.

**Resolution:** these observations drove the decision to drop
`-profile test` entirely — see [task 10](10_ci_stub_only.md).
Fast CI becomes `-stub-run` only; real-tool + biology validation moves
to nightly [`-profile integration`](11_integration_tests.md) with
Tier 2 fixtures. The tiny fixtures described in this task's §5 will
be deleted by task 10; the pre-recruited generation recipe migrates
into [task 11 §6](11_integration_tests.md#6-fixture-generation-recipe)
for the larger integration fixtures.
