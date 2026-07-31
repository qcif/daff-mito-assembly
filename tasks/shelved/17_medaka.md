# Task 17 — P2 stage 8: `MEDAKA`

**Phase:** P2 (from [spec §6](../spec/06-phases.md)).
**Goal:** Replace the P0 stub for stage 8 with a real `medaka_consensus`
invocation that polishes the METAFLYE draft against the gated recruited
FASTQ when `--polish` is passed. This is a thin biocontainer wrapper —
no bespoke Python, no [C1–C7 custom logic](../spec/02-stages.md#22-custom-logic-components).

Per [spec §2 stage 8](../spec/02-stages.md#2-stage-detail) MEDAKA is
**opt-in** (`--polish` boolean flag, default `false`); when disabled the
raw METAFLYE assembly flows straight into `BANDAGE_NG` / `BIN_TARGET`.
Rationale for default-off is spelled out in the same row: Dorado SUP
R10.4.1 reads are already ~Q20+, medaka adds meaningful runtime, and the
cleanup is only useful for marginal-quality ORF calls downstream. The
sample's `metadata.json` / `report.html` record the `polished` flag so
downstream ORF outcomes are attributable ([spec §6a](../spec/06a-reports.md)).

**Prerequisite:** [task 16 — METAFLYE](completed/16_metaflye.md) is real
(it is), so `METAFLYE.out.assembly` carries a real draft
(`assembly.fasta`, `assembly_graph.gfa`, `assembly_info.txt`) into this
stage. The channel wiring already exists in
[main.nf:98-100](../main.nf#L98-L100) —
`params.polish ? MEDAKA(METAFLYE.out.assembly, ch_for_assembly.map { it[1] }).assembly : METAFLYE.out.assembly`
— and requires no topology change, though a small robustness tweak is
suggested in §3 below.

**Exit criteria:**

- `nextflow run . -profile integration --polish true` produces, for
  every `status: "ok"` sample:
    - `results/<sample_id>/assembly/medaka/<sample_id>.polished.fasta`
  A real medaka output, not a `touch`ed placeholder. Total polished
  length is within ±5% of the input METAFLYE assembly length (medaka
  edits homopolymers, does not add/remove contigs).
- MEDAKA runs in the pinned medaka biocontainer
  `quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0` (SHA-pinned in
  [`conf/containers.config`](../conf/containers.config)) — no host tools.
- `-profile stub -stub-run --polish true` still green (stub block
  unchanged; polish branch wired).
- Default integration (`--polish` unset → false) behaviour unchanged:
  MEDAKA is not invoked, downstream still sees `METAFLYE.out.assembly`.
- `bin/` unchanged, `scripts/tests/` unchanged — biocontainer wrapper,
  no C-component logic per [spec §5a rule of thumb](../spec/05-test-data.md#5a-tests).

**Not in scope:**

- Making `--polish` default `true`. Spec §2 stage 8 is explicit: default
  off, opt-in per run. Any future flip requires a spec amendment.
- Model auto-detection. We ship a single Dorado SUP R10.4.1 default; if
  a run's basecaller diverges the operator overrides
  `--medaka_model`. Auto-detection from FASTQ headers is a spec §7 /
  open-question item, not P2.
- Iterative polishing (multiple medaka rounds). Single-pass only —
  matches the CLAW reference workflow and every published organelle
  pipeline surveyed in [spec §7](../spec/07-open-questions.md).
- Polishing the GFA / `assembly_info.txt`. medaka polishes contigs
  only; the graph + per-contig stats pass through unchanged from
  METAFLYE (already the current output shape).

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pinned by SHA, never `latest`.
- No host tools; `medaka_consensus` comes from the container.
- Gated reads + draft assembly enter via `input:` — no bare
  `${params.<file>}`.

---

## 1. Params — `nextflow.config`

Two additions alongside the existing `polish` boolean at
[nextflow.config:45](../nextflow.config#L45):

```groovy
// --- Assembly ---
polish        = false

// Medaka model — matched to Dorado SUP R10.4.1 (spec §2 stage 8).
// Override at runtime if the basecaller / model differs.
medaka_model  = 'r1041_e82_400bps_sup_v5.0.0'
```

**Model selection.** Confirm the exact model string against
`medaka tools list_models` inside the pinned container before locking
in — medaka's model catalog rotates faster than the tool version, and
2.2.2 may ship a newer default. Record the resolved model name in
Outcomes. The convention `r1041_e82_400bps_sup_v<X>` tracks Dorado's
R10.4.1 SUP basecaller family; pick the newest `v*` that 2.2.2 knows
about.

## 2. Container pin — `conf/containers.config`

Replace the current `python:3.12-slim` stub + `TODO P2` comment for
`MEDAKA` (at
[conf/containers.config:52-55](../conf/containers.config#L52-L55)) with
the pinned medaka biocontainer:

```groovy
withName: 'MEDAKA' {
    container = 'quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0@sha256:<pinned>'
}
```

**Resolving the SHA:** pull the tag locally (mind the CLAUDE.md
`--user` habit for interactive runs), then read the digest:

```bash
docker pull quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0
docker inspect --format='{{index .RepoDigests 0}}' \
    quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0
```

Record the exact resolved digest in Outcomes. Single-tool biocontainer
— no `mulled-build` needed
([spec §1a container selection order](../spec/01-pipeline-flow.md#1a-engineering-constraints),
step 1).

## 3. `modules/local/medaka.nf` (replace stub)

Replace the current stub `script:` block (at
[modules/local/medaka.nf:25-31](../modules/local/medaka.nf#L25-L31))
with a real `medaka_consensus` invocation. Retain the `stub:` block
unchanged for `-profile stub`.

Pseudocode of the intended `script:` block:

```groovy
script:
def model = params.medaka_model
"""
medaka_consensus \\
    -i ${reads} \\
    -d ${assembly} \\
    -o medaka_out \\
    -m ${model} \\
    -t ${task.cpus}
mv medaka_out/consensus.fasta ${meta.sample_id}.polished.fasta
"""
```

Notes:

- **Output shape unchanged.** The module already emits
  `tuple(meta, ${meta.sample_id}.polished.fasta, gfa, info)` — matches
  the downstream `ch_assembly` fork at
  [main.nf:98-100](../main.nf#L98-L100). GFA + info pass through from
  the input tuple; medaka does not touch them. No `output:` change.
- **`-o medaka_out` then move.** Same rationale as METAFLYE
  ([task 16 §3](completed/16_metaflye.md)): flatten the one file we
  need into the task work-dir. Everything else medaka writes (BAM,
  HDF5 intermediates, logs) is discarded — reproducibility is via the
  pinned container + params.
- **`process_high` label retained.** Medaka is GPU-optional but CPU
  works within the CI ceiling (2 vCPU / 14 GB / 60 m on
  ubuntu-latest) for the three Tier 2 organelle-scale fixtures.
- **Channel-pairing note.** The current wiring at
  [main.nf:99](../main.nf#L99) pairs
  `METAFLYE.out.assembly` with `ch_for_assembly.map { it[1] }` by
  positional cardinality — fine for a single-fixture run but fragile
  under reordering. If a per-meta join is trivial (both channels are
  keyed on `meta`), tighten it to `join(by: 0)` while touching the
  wiring; if not, leave a comment noting the ordering assumption.

## 4. Integration-test wiring

Two considerations, one intentional decision.

**Default integration run (`--polish` unset).** No change. MEDAKA is
not invoked; existing 14/14 assertions
([task 16 outcomes](completed/16_metaflye.md#10-outcomes)) still pass
unchanged.

**Opt-in polish smoke test.** Add a targeted assertion block in
[`tests/integration/assertions.sh`](../tests/integration/assertions.sh),
guarded on the presence of the polished FASTA (so it runs only when
the CI invocation passes `--polish`, and stays silent otherwise). This
avoids branching the integration workflow while catching medaka
regressions when we choose to exercise them.

Sketch:

```bash
# MEDAKA polish (only if --polish was passed to the integration run)
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    polished="$OUTDIR/$sample/assembly/medaka/${sample}.polished.fasta"
    draft="$OUTDIR/$sample/assembly/assembly.fasta"
    [[ -s "$polished" ]] || continue   # not asked for; skip silently
    p_bp=$(grep -v '^>' "$polished" | tr -d '\n' | wc -c)
    d_bp=$(grep -v '^>' "$draft"    | tr -d '\n' | wc -c)
    # Medaka edits homopolymers; total length within ±5% of draft.
    delta_pct=$(( (p_bp - d_bp) * 100 / d_bp ))
    if (( delta_pct < -5 || delta_pct > 5 )); then
        echo "FAIL: $sample polish delta ${delta_pct}% outside ±5%"
        FAILED=1
    else
        echo "OK:   $sample polish delta ${delta_pct}% (${p_bp} bp vs ${d_bp} bp draft)"
    fi
done
```

Update the progressive-uncomment plan header at the top of
`assertions.sh` to add:

```
#   MEDAKA (real polisher)        | Polished length within ±5% of draft (guarded on presence) [done]
```

**CI invocation.** The scheduled `Integration` workflow keeps its
default (no `--polish`). A separate manual verification during task
execution runs `nextflow run . -profile integration --polish true`
locally to prove the wiring and record delta numbers in Outcomes; no
CI matrix expansion needed. If polish becomes a supported production
knob, add an explicit `Integration (polish)` workflow job in a
follow-up task.

**No new fixtures.** `assembly.fasta` from METAFLYE is the input; the
delta assertion is derived from the same file. No `expected/*/*.json`
changes.

## 5. Fast CI

No change. The stub `script:` block still `touch`es the polished FASTA;
Nextflow container-coverage check is satisfied by the SHA-pinned tag;
`-stub-run` with or without `--polish` continues to pass without
pulling the medaka image.

Verify: `nextflow run . -profile stub -stub-run --polish true` — 18/18
processes hit their stub (17 previously + MEDAKA now that the polish
branch fires), `tests/output/STUB-01/` unchanged.

## 6. Unit tests

**None required.** Biocontainer wrapper with no
[C-component custom logic](../spec/02-stages.md#22-custom-logic-components).
Per [spec §5a](../spec/05-test-data.md#5a-tests):

- Channel wiring covered by `-stub-run` (fast CI on every push,
  extended to include `--polish true`).
- Command-line correctness + biology covered by manual
  `-profile integration --polish true` (this task) plus the guarded
  assertion block for any future scheduled polish run.

If future work promotes model selection into a helper
(e.g. `bin/pick_medaka_model.py` reading FASTQ headers), a
`scripts/tests/` module gets added then — not now.

## 7. Verification

1. `nextflow run . -profile stub -stub-run` — 17/17 stubs (polish off,
   MEDAKA not invoked). Unchanged from task 16.
2. `nextflow run . -profile stub -stub-run --polish true` — 18/18
   stubs; MEDAKA fires and its stub `touch`es
   `${sample_id}.polished.fasta`.
3. `nextflow run . -profile integration` locally — assembly outputs +
   14 assertions unchanged from task 16.
4. `nextflow run . -profile integration --polish true` locally
   (real medaka on all three Tier 2 fixtures):
   - `INT-ANIMAL-01/assembly/medaka/INT-ANIMAL-01.polished.fasta` —
     ~349 kb ±5% (matches
     [task 16 outcomes](completed/16_metaflye.md#10-outcomes) input
     size).
   - `INT-PLANT-01-pt/assembly/medaka/INT-PLANT-01-pt.polished.fasta`
     — ~155 kb ±5%.
   - `INT-PLANT-01-mt/assembly/medaka/INT-PLANT-01-mt.polished.fasta`
     — ~270 kb ±5%.
5. `bash tests/integration/assertions.sh` — every pre-existing
   assertion green; guarded MEDAKA block reports OK for all three
   samples when polish was on, silent otherwise.
6. `flake8` — no Python touched; nothing to run.
7. CI green on `Tests` (fast, unchanged) + nightly `Integration`
   (unchanged, polish off) on the PR.

## 8. Deliverables checklist

- [ ] [`nextflow.config`](../nextflow.config) — `params.medaka_model`
      added; existing `params.polish = false` unchanged.
- [ ] [`conf/containers.config`](../conf/containers.config) — SHA-pinned
      medaka biocontainer for `MEDAKA`; `TODO P2` comment removed.
- [ ] [`modules/local/medaka.nf`](../modules/local/medaka.nf) — real
      `script:` block invoking `medaka_consensus`; stub retained.
- [ ] Optional: tighten
      [main.nf:99](../main.nf#L99) channel pairing to `join(by: 0)`
      *or* add a comment noting the ordering assumption.
- [ ] [`tests/integration/assertions.sh`](../tests/integration/assertions.sh)
      — guarded polish-delta block added; progressive-uncomment plan
      header updated.
- [ ] Fast CI (`Tests`) + `Integration` workflows both green on the PR
      (polish default unchanged; wiring proven via local
      `--polish true` run recorded in Outcomes).

## 9. Notes / non-issues

- **No channel-topology change.** Fork at
  [main.nf:98-100](../main.nf#L98-L100) already handles both
  polish/no-polish; downstream `ch_assembly` is agnostic to the
  branch. This task swaps one stub for one real invocation.
- **`polished` provenance.** COLLATE eventually needs to record the
  effective polish state in `metadata.json` per
  [spec §2 stage 8](../spec/02-stages.md#2-stage-detail) /
  [spec §6a](../spec/06a-reports.md). That's a COLLATE task
  concern (task TBD in P4), not MEDAKA's job here — MEDAKA just
  produces the polished FASTA (or doesn't fire).
- **Model version churn.** Medaka's model catalog changes independent
  of the pip release. If `medaka tools list_models` inside the
  pinned image no longer lists `r1041_e82_400bps_sup_v5.0.0`, pick the
  newest R10.4.1 SUP model it does list and record in Outcomes.
- **Resource ceiling.** Medaka fits comfortably within the 2 vCPU /
  14 GB / 60 m ceiling for organelle-scale drafts (largest fixture is
  ~349 kb from `--meta`). If a `plant_mt` fragmentation regression
  balloons the draft, medaka slows proportionally — address only if
  it actually OOMs / times out.
- **GPU not required.** Medaka 2.x runs its neural net on CPU by
  default; the biocontainer is CPU-only. Do not attempt GPU wiring in
  this task.

## 10. Outcomes

_To be filled in on completion._
