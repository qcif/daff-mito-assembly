# Task 16 — P2 stage 7: `METAFLYE`

**Phase:** P2 (from [spec §6](../spec/06-phases.md)).
**Goal:** Replace the P0 stub for stage 7 with a real Flye `--meta`
invocation that assembles the gated recruited FASTQ into a draft
organelle assembly, an assembly graph, and Flye's per-contig
`assembly_info.txt`. This is a thin biocontainer wrapper — no bespoke
Python, no [C1–C7 custom logic](../spec/02-stages.md#22-custom-logic-components).

Per [spec §2 stage 7](../spec/02-stages.md#2-stage-detail) Flye is invoked
in metagenomic mode (`--meta`) with the ONT-HQ read model (`--nano-hq`),
a per-target `--genome-size` hint ([spec §3.4](../spec/03-organelles.md#34-flye---genome-size-hint-per-target)),
The three outputs (`assembly.fasta`, `assembly_graph.gfa`,
`assembly_info.txt`) all feed [BIN_TARGET (task TBD)](../spec/02-stages.md#22-custom-logic-components);
`assembly_info.txt` in particular is **load-bearing** for the
coverage-spike component of BIN_TARGET and must not be dropped.

**Prerequisite:** [task 15 — COVERAGE_GATE](completed/15_coverage_gate.md)
is real (it is), so `COVERAGE_GATE.out.gated`'s `ok` branch carries a
real, coverage-capped gzipped FASTQ into this stage. The channel wiring
already exists in [main.nf:96](../main.nf#L96) — `METAFLYE(ch_for_assembly)`
— and requires no topology change.

**Exit criteria:**

- `-profile integration` produces, for every `status: "ok"` sample:
    - `results/<sample_id>/assembly/assembly.fasta`
    - `results/<sample_id>/assembly/assembly_graph.gfa`
    - `results/<sample_id>/assembly/assembly_info.txt`
  All three are real Flye outputs, not `touch`ed placeholders. The
  FASTA contains at least one contig; the GFA has ≥ 1 `S` line;
  `assembly_info.txt` has a header + one row per contig.
- Assembly-length assertion block in
  [`tests/integration/assertions.sh:71-96`](../tests/integration/assertions.sh#L71-L96)
  uncommented and green: total assembly length for each Tier 2 fixture
  falls inside the per-target bounds recorded in
  [`tests/integration/expected/<target>/assembly_bounds.json`](../tests/integration/expected/).
- METAFLYE runs in the pinned Flye biocontainer
  `quay.io/biocontainers/flye:2.9.6--py311h93bbee8_1` (SHA-pinned in
  [`conf/containers.config`](../conf/containers.config)) — no host tools.
- `-profile stub -stub-run` still green (stub block unchanged).
- `bin/` unchanged, `scripts/tests/` unchanged — this is an
  off-the-shelf tool wrapper per [spec §5a rule of thumb](../spec/05-test-data.md#5a-tests);
  channel wiring is covered by `-stub-run`, biology by nightly
  integration.

**Not in scope:**

- `MEDAKA` polish stage (P2 sibling task) — remains stub. The
  `ch_assembly = params.polish ? MEDAKA(...) : METAFLYE.out.assembly`
  branch at [main.nf:98](../main.nf#L98) already handles both.
- `BANDAGE_NG` real-tool wiring (P2 sibling task).
- `--nano-corr` mode selection ([spec §9 item 5](../spec/07-open-questions.md)) —
  default `--nano-hq` per spec until an error-correction stage is
  introduced upstream. Revisit when that decision lands.
- `plant_mt` `--genome-size` re-evaluation ([spec §3.4](../spec/03-organelles.md#34-flye---genome-size-hint-per-target)
  "Re-evaluate the `plant_mt` hint in P2 if mt assembly fragments") —
  ship the `2m` default; if the integration fixture fragments badly,
  file a follow-up task in [tasks/todo.md](todo.md).
- Circularity check on `animal_mt` — that's a BIN_TARGET concern
  ([spec §3.3](../spec/03-organelles.md#33-specific-issues-and-decisions)),
  not METAFLYE's.

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pinned by SHA, never `latest`.
- No host tools; `flye` comes from the container.
- Gated reads enter via `input: path` — no bare `${params.<file>}`.

---

## 1. Per-target Flye parameters — `nextflow.config`

Add an `assembly_size_hints` table alongside `coverage_limits` in
[`nextflow.config`](../nextflow.config). Kept separate from
`coverage_limits.nominal_size` because the two serve different
purposes: `nominal_size` is the coverage-estimation denominator
(exact-ish per-target constant), while Flye's `--genome-size` is a
loose hint tuned for the assembler's coverage estimator
([spec §3.4](../spec/03-organelles.md#34-flye---genome-size-hint-per-target)).

```groovy
// Flye --genome-size hints — spec §3.4.
// Distinct from coverage_limits.nominal_size: this is an advisory hint
// to Flye's coverage estimator, not a length filter.
assembly_size_hints = [
    animal_mt: '20k',
    plant_pt:  '150k',
    plant_mt:  '2m',
]
```

## 2. Container pin — `conf/containers.config`

Replace the current `python:3.12-slim` stub + `TODO P2` comment for
`METAFLYE` (at [conf/containers.config:49-52](../conf/containers.config#L49-L52))
with the pinned Flye biocontainer:

```groovy
withName: 'METAFLYE' {
    container = 'quay.io/biocontainers/flye:2.9.6--py311h93bbee8_1@sha256:<pinned>'
}
```

**Resolving the SHA:** pull the tag locally, then read the digest:

```bash
docker pull quay.io/biocontainers/flye:2.9.6--py311h93bbee8_1
docker inspect --format='{{index .RepoDigests 0}}' \
    quay.io/biocontainers/flye:2.9.6--py311h93bbee8_1
```

Record the exact resolved digest in the Outcomes section on completion.
Single-tool biocontainer — no `mulled-build` needed
([spec §1a container selection order](../spec/01-pipeline-flow.md#1a-engineering-constraints),
step 1).

## 3. `modules/local/metaflye.nf` (replace stub)

Replace the current stub script block (at
[modules/local/metaflye.nf:24-29](../modules/local/metaflye.nf#L24-L29))
with a real Flye invocation. Retain the stub block unchanged for
`-profile stub`.

Pseudocode of the intended `script:` block:

```groovy
script:
def genome_size = params.assembly_size_hints[meta.assembly_target]
"""
flye --meta --nano-hq ${reads} \\
     --genome-size ${genome_size} \\
     --threads ${task.cpus} \\
     --out-dir flye_out
mv flye_out/assembly.fasta       assembly.fasta
mv flye_out/assembly_graph.gfa   assembly_graph.gfa
mv flye_out/assembly_info.txt    assembly_info.txt
"""
```

Notes:

- **Output shape unchanged.** The module already emits
  `tuple(meta, assembly.fasta, assembly_graph.gfa, assembly_info.txt)`
  — matches downstream expectations at
  [main.nf:98-102](../main.nf#L98-L102). No `output:` change.
- **`--out-dir` then move.** Flye insists on writing to a directory it
  creates; we flatten the three files we care about into the task
  work-dir so `output: path("assembly*")` matches without a subdir.
  Everything else Flye writes (params.json, log, intermediate graphs)
  is discarded — reproducibility is via the pinned container + params,
  not Flye's per-run bookkeeping.
- **`process_high` label retained.** Flye is CPU + memory heavy; the
  `resourceLimits` in [`conf/integration.config`](../conf/integration.config)
  (2 vCPU / 14 GB / 60 m on ubuntu-latest) is what the CI runner can
  provide. Skim-depth organelle assembly fits inside this budget for
  the three Tier 2 fixtures ([spec §5](../spec/05-test-data.md)) —
  confirm on first real run.
- **Read mode `--nano-hq`.** Assumes Dorado SUP R10.4.1 input; if an
  error-correction stage lands upstream, switch to `--nano-corr`
  ([spec §9 item 5](../spec/07-open-questions.md)). Out of scope here.

## 4. Integration-test wiring

Uncomment the METAFLYE assertion block in
[`tests/integration/assertions.sh:71-96`](../tests/integration/assertions.sh#L71-L96),
adjusting only the assembly path if needed to match the publish
directory (`<outdir>/<sample_id>/assembly/assembly.fasta`, not
`<outdir>/<sample_id>/organelle_assembly.fasta` — that's a
COLLATE-emitted path). Also remove the TODO comment on the
`METAFLYE (real assembler)` line in the progressive-uncomment plan
header at the top of the file.

Sketch of the block after adjustment:

```bash
declare -A SAMPLE_TARGET=(
    [INT-ANIMAL-01]=animal_mt
    [INT-PLANT-01-pt]=plant_pt
    [INT-PLANT-01-mt]=plant_mt
)
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    asm="$OUTDIR/$sample/assembly/assembly.fasta"
    if [[ ! -s "$asm" ]]; then
        echo "FAIL: $sample assembly missing or empty"
        FAILED=1
        continue
    fi
    asm_bp=$(grep -v '^>' "$asm" | tr -d '\n' | wc -c)
    target="${SAMPLE_TARGET[$sample]}"
    bounds="tests/integration/expected/${target}/assembly_bounds.json"
    min_bp=$(jq .min_bp "$bounds")
    max_bp=$(jq .max_bp "$bounds")
    if (( asm_bp < min_bp || asm_bp > max_bp )); then
        echo "FAIL: $sample assembly ${asm_bp} bp outside [${min_bp}, ${max_bp}]"
        FAILED=1
    else
        echo "OK:   $sample assembly ${asm_bp} bp in [${min_bp}, ${max_bp}]"
    fi
done
```

**Progressive-assertion note.** Only the total-length bounds are
checked — no per-contig count, no circularity, no annotation content.
Those belong to BIN_TARGET / MITOS2 / MINIPROT_EXTRACT assertion blocks
and stay commented until their stages land. `assembly_bounds.json`
already exists for all three targets (see
[`tests/integration/expected/*/assembly_bounds.json`](../tests/integration/expected/)) —
no new fixture generation.

**Publish-path scoping.** The module currently publishes under
`<outdir>/<sample_id>/assembly/` gated on
`params.publish_intermediates`. Integration runs need those three files
available for the assertion, so either:

1. Enable `publish_intermediates = true` in
   [`conf/integration.config`](../conf/integration.config) (broadest —
   also publishes recruit / gated FASTQs; probably desirable for
   integration debugging anyway); or
2. Add a second unconditional `publishDir` for
   `pattern: 'assembly.fasta'` (or all three files) — mirrors the
   pattern used by [COVERAGE_GATE for
   `sample_status.json`](../modules/local/coverage_gate.nf#L12-L13).

Pick (1) if no other integration outputs currently need special-cased
publish; (2) if we want to keep the intermediate flag opt-in. Record
the choice in Outcomes.

## 5. Fast CI

No change. The stub script block still `touch`es the three files, the
Nextflow container-coverage check remains satisfied (SHA-pinned tag,
no `latest`), and `-stub-run` continues to pass without pulling the
Flye image.

## 6. Unit tests

**None required.** This stage is a biocontainer wrapper with no
[C-component custom logic](../spec/02-stages.md#22-custom-logic-components).
Per [spec §5a](../spec/05-test-data.md#5a-tests):

- Channel wiring covered by `-stub-run` (fast CI on every push).
- Command-line correctness + biology covered by nightly
  `-profile integration` + `assertions.sh`.

If future work promotes any of the per-target parameter selection into
a helper (e.g. a `bin/flye_args.py`), a `scripts/tests/` module gets
added then — not now.

## 7. Verification

1. `nextflow run . -profile stub -stub-run` — 17/17 processes hit
   their stub, `tests/output/STUB-01/` unchanged.
2. `nextflow run . -profile integration` locally (needs the fetched
   fixtures + refdata bundle from
   [task 14](completed/14_complete_integration_ci.md)):
   - `INT-ANIMAL-01/assembly/assembly.fasta` — 1 contig, total
     14–20 kb (matches
     [`expected/animal_mt/assembly_bounds.json`](../tests/integration/expected/animal_mt/assembly_bounds.json)).
   - `INT-PLANT-01-pt/assembly/assembly.fasta` — 1–4 contigs, total
     145–165 kb (matches
     [`expected/plant_pt/assembly_bounds.json`](../tests/integration/expected/plant_pt/assembly_bounds.json)).
   - `INT-PLANT-01-mt/assembly/assembly.fasta` — plant mt is
     multi-contig; total 200–700 kb (matches
     [`expected/plant_mt/assembly_bounds.json`](../tests/integration/expected/plant_mt/assembly_bounds.json)).
   - Every `assembly_info.txt` has a header + ≥ 1 data row.
   - Every `assembly_graph.gfa` has ≥ 1 `S` line.
3. `bash tests/integration/assertions.sh` — all pre-existing
   assertions still green + the newly uncommented METAFLYE block
   reports OK for all three samples.
4. `flake8` — no Python touched; nothing to run.
5. CI green on `Tests` (fast) + nightly `Integration` on the PR.

## 8. Deliverables checklist

- [x] [`nextflow.config`](../nextflow.config) — `params.assembly_size_hints`
      table added.
- [x] [`conf/containers.config`](../conf/containers.config) — SHA-pinned
      Flye biocontainer for `METAFLYE`, `TODO P2` comment removed.
- [x] [`modules/local/metaflye.nf`](../modules/local/metaflye.nf) —
      real `script:` block invoking Flye per-target; stub retained.
- [x] Publish path for `assembly.fasta` visible under
      `<outdir>/<sample_id>/assembly/` during
      `-profile integration` runs (enabled `publish_intermediates = true`
      in `conf/integration.config`).
- [x] [`tests/integration/assertions.sh`](../tests/integration/assertions.sh)
      — METAFLYE assembly-length block uncommented; header TODO line
      for METAFLYE removed.
- [ ] Fast CI (`Tests`) + `Integration` workflows both green on the PR.

## 9. Notes / non-issues

- **No channel-topology change.** [main.nf:96](../main.nf#L96) already
  wires `METAFLYE(ch_for_assembly)`; downstream `ch_assembly` handles
  the polish/no-polish fork; BANDAGE_NG / BIN_TARGET already consume
  METAFLYE outputs (all still stubs). This task swaps one stub for
  one real invocation — nothing else moves.
- **No new fixtures.** `assembly_bounds.json` for all three targets
  was created in [task 11](completed/11_integration_tests.md) /
  [task 12](completed/12_azure_integration_fixtures.md). If a fixture
  regenerates during work here, update the bounds file rather than
  loosening the assertion.
- **`--meta` retained.** Even though each sample-sheet row is one
  organelle, `--meta` tolerates low-level contamination
  ([brief.md §3.5](../brief.md)) — same rationale as the spec entry
  for stage 7.
- **`plant_mt` fragmentation risk.** If the `plant_mt` fixture
  assembles into many small pieces below the 200 kb lower bound, that
  is a spec-flagged known concern
  ([spec §3.4](../spec/03-organelles.md#34-flye---genome-size-hint-per-target),
  [§3.3](../spec/03-organelles.md#33-specific-issues-and-decisions)).
  Not a task blocker — file a P2 follow-up to sweep `--genome-size` /
  `--min-overlap`, and consider loosening the plant_mt lower bound
  in `assembly_bounds.json` if the current 200 kb is too tight for
  the SRR11315861 fixture at skim depth.
- **Resource ceiling.** If Flye OOMs under the 14 GB CI ceiling on
  `plant_mt`, options: (a) reduce `plant_mt` fixture size (COVERAGE_GATE
  will already cap it at 300× — verify the effective input is < a few
  hundred MB); (b) split integration into per-sample matrix jobs.
  Address only if it actually OOMs — spec [§5](../spec/05-test-data.md)
  sizes these fixtures deliberately small.

## 10. Outcomes

- Container: `quay.io/biocontainers/flye:2.9.6--py311h93bbee8_1@sha256:eb57665b9d12c43f1112cfd7af64ffd94236a815ff1794298088b01500632b9f`
- **`--asm-coverage` dropped**: incompatible with `--meta` in Flye 2.9.6.
  Spec §3.5 removed accordingly. `--meta` is required to tolerate
  contamination and takes precedence.
- **`publish_intermediates = true` in integration.config**: simplest
  way to expose `assembly/assembly.fasta` to the assertion script; also
  useful for integration debugging of all upstream stages.
- **`animal_mt` assembly bounds widened** (10 kb – 2 Mb): `--meta` mode
  assembles off-target survivors alongside the mitochondrion (15 contigs,
  ~349 kb total). The circular 71× contig_8 (~17 kb) is the true mt;
  BIN_TARGET will select it. The METAFLYE assertion is a coarse sanity
  check; the tight biology check moves to a future BIN_TARGET assertion.
- Integration results: INT-ANIMAL-01 ~349 kb (15 contigs, 2 circular at
  ~17 kb); INT-PLANT-01-pt 155 kb; INT-PLANT-01-mt 270 kb.
- All 14 integration assertions pass; 17/17 stub-run green.
