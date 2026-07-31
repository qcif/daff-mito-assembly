# Task 17 — P2 stage 9: `BANDAGE_NG`

**Phase:** P2 (from [spec §6](../spec/06-phases.md)).
**Goal:** Replace the P0 stub for stage 9 with a real
[BandageNG](https://github.com/asl/BandageNG) image render of the
METAFLYE assembly graph. Diagnostic-only stage — the PNG is a
human-review artefact that lands in the per-sample bundle. This is a
thin biocontainer wrapper — no bespoke Python, no
[C1–C7 custom logic](../spec/02-stages.md#22-custom-logic-components).

Per [spec §2 stage 9](../spec/02-stages.md#2-stage-detail) BandageNG
consumes `assembly_graph.gfa` and emits a PNG of the graph, which is
consumed downstream by `COLLATE` for inclusion in the per-sample
`report.html` (see [spec §6a.2](../spec/06a-reports.md)).

**Prerequisite:** [task 16 — METAFLYE](completed/16_metaflye.md) is real,
so `METAFLYE.out.assembly` carries a real `assembly_graph.gfa` (with ≥ 1
`S` line) into this stage. Channel wiring already exists in
[main.nf:102](../main.nf#L102) — `BANDAGE_NG(ch_assembly)` — and requires
no topology change.

**Exit criteria:**

- `-profile integration` produces, for every `status: "ok"` sample:
    - `results/<sample_id>/assembly/<sample_id>.graph.png` — a real
      PNG rendered by BandageNG (not a `touch`ed placeholder). File
      size > 1 kB and identifiable as PNG by magic bytes.
- BANDAGE_NG runs in the pinned BandageNG biocontainer
  `quay.io/biocontainers/bandage_ng:<version>` (SHA-pinned in
  [`conf/containers.config`](../conf/containers.config#L56-L59)) — no
  host tools.
- Integration assertion block added to
  [`tests/integration/assertions.sh`](../tests/integration/assertions.sh)
  checking PNG magic bytes for each sample.
- `-profile stub -stub-run` still green (stub block unchanged).
- `bin/` unchanged, `scripts/tests/` unchanged — off-the-shelf tool
  wrapper per [spec §5a rule of thumb](../spec/05-test-data.md#5a-tests).

**Not in scope:**

- Layout or styling tuning of the graph (edge colouring by depth,
  node labels). BandageNG defaults are sufficient for diagnostic use.
- SVG output. Report renderer (P4) can consume PNG directly; adding
  SVG doubles publish surface for no downstream consumer today.
- Interactive HTML embed of the graph. `report.html` links or embeds
  the PNG per [spec §6a.2](../spec/06a-reports.md).
- Any reaction to graph *content* (edge count, canonicalisation) —
  that is `BIN_TARGET`'s job ([task 18](18_bin_target.md),
  [spec §3.6](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation)).

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pinned by SHA, never `latest`.
- No host tools; `Bandage` binary comes from the container.
- GFA enters via `input: path` — no bare `${params.<file>}`.

---

## 1. Container pin — `conf/containers.config`

Replace the `python:3.12-slim` stub + `TODO P2` comment for
`BANDAGE_NG` (at [conf/containers.config:56-59](../conf/containers.config#L56-L59))
with the pinned BandageNG biocontainer already staged in the TODO
comment:

```groovy
withName: 'BANDAGE_NG' {
    container = 'quay.io/biocontainers/bandage_ng:2026.6.1--hca0ed12_0@sha256:<pinned>'
}
```

**Resolving the SHA:** as in
[task 16 §2](completed/16_metaflye.md), pull the tag locally then read
the digest:

```bash
docker pull quay.io/biocontainers/bandage_ng:2026.6.1--hca0ed12_0
docker inspect --format='{{index .RepoDigests 0}}' \
    quay.io/biocontainers/bandage_ng:2026.6.1--hca0ed12_0
```

Record the resolved digest in Outcomes on completion. Single-tool
biocontainer — no `mulled-build` needed
([spec §1a container selection order](../spec/01-pipeline-flow.md#1a-engineering-constraints)).

**Note on `latest`.** If Bioconda's newest tag has shifted since this
task was drafted, pin whatever the current pinned build is — the
constitution's [rule 11](../CONSTITUTION.md) forbids `latest` but does
not mandate a specific version. Record the chosen tag in Outcomes.

## 2. `modules/local/bandage_ng.nf` (replace stub)

Replace the current stub script block (at
[modules/local/bandage_ng.nf:22-27](../modules/local/bandage_ng.nf#L22-L27))
with a real `Bandage image` invocation. Retain the stub block unchanged
for `-profile stub`.

Pseudocode of the intended `script:` block:

```groovy
script:
"""
Bandage image ${gfa} ${meta.sample_id}.graph.png \\
    --height 600 \\
    --width 800
"""
```

Notes:

- **Output shape unchanged.** The module already emits
  `tuple(meta, assembly, gfa, info, graph_png)` — matches
  `BIN_TARGET`'s expected input at
  [main.nf:105](../main.nf#L105). No `output:` change; downstream
  wiring is unaffected.
- **Fixed dimensions.** 800×600 keeps PNG size predictable for the
  per-sample report layout; BandageNG scales the graph to fit.
- **Headless.** The `Bandage image` subcommand does not need an X
  display; the biocontainer ships without xvfb.
- **Process label unchanged.** `process_low` is correct — graph
  rendering on organelle-scale GFAs is CPU-cheap and memory-cheap.

## 3. Integration-test wiring

Add a BANDAGE_NG assertion block after the METAFLYE block in
[`tests/integration/assertions.sh`](../tests/integration/assertions.sh):

```bash
# BANDAGE_NG is real (task 17):
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    png="$OUTDIR/$sample/assembly/${sample}.graph.png"
    if [[ ! -s "$png" ]]; then
        echo "FAIL: $sample graph PNG missing or empty"
        FAILED=1
        continue
    fi
    # PNG magic bytes: 89 50 4E 47 0D 0A 1A 0A
    magic=$(head -c 8 "$png" | xxd -p)
    if [[ "$magic" != "89504e470d0a1a0a" ]]; then
        echo "FAIL: $sample graph PNG magic bytes bad ($magic)"
        FAILED=1
    else
        echo "OK:   $sample graph PNG ($(wc -c < "$png") bytes)"
    fi
done
```

Also update the progressive-uncomment header at the top of the file to
tick off the BANDAGE_NG row.

**No new fixtures.** The PNG is asserted structurally (magic bytes,
non-empty) — the graph content itself depends on assembly output, which
already has its own bounds in
[`assembly_bounds.json`](../tests/integration/expected/). No
`bandage_bounds.json` is needed; the render is diagnostic, not
biology-checkable in a stable way.

**Publish path.** The module publishes to
`<outdir>/<sample_id>/assembly/` gated on `params.publish_intermediates`.
That flag is already `true` in
[`conf/integration.config`](../conf/integration.config#L12) since
task 16, so the PNG will land where the assertion looks for it.

## 4. Fast CI

No change. The stub script block still `touch`es the PNG, the
container-coverage lint stays satisfied (SHA-pinned tag, no `latest`),
and `-stub-run` continues to pass without pulling the BandageNG image.

## 5. Unit tests

**None required.** Off-the-shelf tool wrapper with no
[C-component custom logic](../spec/02-stages.md#22-custom-logic-components).
Per [spec §5a](../spec/05-test-data.md#5a-tests):

- Channel wiring covered by `-stub-run` (fast CI on every push).
- Real-tool behaviour covered by nightly `-profile integration` +
  `assertions.sh`.

## 6. Verification

1. `nextflow run . -profile stub -stub-run` — 17/17 processes hit
   their stub, no regressions in `tests/output/STUB-01/`.
2. `nextflow run . -profile integration` locally:
   - Each `<sample_id>.graph.png` is a valid PNG > 1 kB.
   - `file` command reports "PNG image data".
   - Manual eyeball: pea aphid PNG shows a small circular graph
     (~15 nodes); plant plastid PNG shows the classic three-node
     quadripartite structure (LSC/IR/SSC) when the IR is resolved.
3. `bash tests/integration/assertions.sh` — new BANDAGE_NG block
   reports OK for all three samples.
4. `flake8` — no Python touched; nothing to run.
5. CI green on `Tests` (fast) + nightly `Integration` on the PR.

## 7. Deliverables checklist

- [x] [`conf/containers.config`](../conf/containers.config) — SHA-pinned
      BandageNG biocontainer for `BANDAGE_NG`, `TODO P2` comment
      removed.
- [x] [`modules/local/bandage_ng.nf`](../modules/local/bandage_ng.nf) —
      real `BandageNG image` `script:` block; stub retained.
- [x] [`tests/integration/assertions.sh`](../tests/integration/assertions.sh) —
      BANDAGE_NG PNG-magic assertion block added; progressive-uncomment
      header updated.
- [ ] Fast CI (`Tests`) + `Integration` workflows both green on the PR.

## 8. Notes / non-issues

- **No channel-topology change.** [main.nf:102](../main.nf#L102)
  already wires `BANDAGE_NG(ch_assembly)`; the emit tuple shape is
  identical to the stub's; downstream `BIN_TARGET` still receives
  `(meta, assembly, gfa, info, graph_png)`.
- **`plant_mt` graph size.** Plant mitogenomes can produce large,
  tangled multipartite graphs
  ([spec §3.3](../spec/03-organelles.md#33-specific-issues-and-decisions)).
  If the PNG becomes unreadable on the SRR11315861 fixture, tune
  `--height` / `--width` up rather than filtering nodes — the graph
  is diagnostic and complexity is the signal.
- **`--iterations` / layout stability.** BandageNG uses a force-directed
  layout; the same GFA will render slightly differently across runs.
  This is fine — the assertion is on file structure, not pixel content.
- After running integration test, prompt the user to review the output PNG/SVG
  manually.

## 9. Outcomes

- Container: `quay.io/biocontainers/bandage_ng:2026.6.1--hca0ed12_0@sha256:dfc461d4795bc7186a6549ec1fe4dfaf23c861c582add29b70a4bae34dbaf1fc`
- Chosen BandageNG version: `2026.6.1`
- **Executable name deviation:** the biocontainer ships the binary as `BandageNG`
  (not `Bandage` as written in the task pseudocode). Updated `script:` block and
  task pseudocode accordingly.
- `-stub-run` green: 17/17 processes.
