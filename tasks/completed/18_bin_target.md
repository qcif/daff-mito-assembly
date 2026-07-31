# Task 18 — P3 stage 10: `BIN_TARGET`

**Phase:** P3 (from [spec §6](../spec/06-phases.md)).
**Goal:** Replace the P0 stub for stage 10 with a real per-contig
binning stage driven by [`bin/bin_target.py`](../bin/) (custom-logic
component C3 in
[spec §2.2](../spec/02-stages.md#22-custom-logic-components)). The
stage classifies each METAFLYE contig as target vs off-target using
the **intersection** of three signals:

1. Coverage spike vs. background (from `assembly_info.txt`).
2. Reference identity to `${assembly_target}.mmi` (minimap2).
3. ORF integrity under the target-appropriate NCBI genetic-code table.

The dominant target contig is emitted as `target.fasta`; the rest are
logged in `secondaries.tsv`. Per
[spec §3.3](../spec/03-organelles.md#33-specific-issues-and-decisions)
the `animal_mt` branch additionally runs an **end-overlap circularity
check** on the selected contig and records the result in the
diagnostics.

**C4 (`plastid_canonicalise.py`) is broken out into two sibling
tasks under a clean-room-firewall design** (upstream ptGAUL is
unlicensed):

- [19_plastid_algorithm_spec.md](19_plastid_algorithm_spec.md) —
  read `combine_gfa.py` + spec §3.6, produce
  `spec/plastid-canonicalisation.md` (algorithm design doc,
  no runnable code).
- [20_plastid_canonicalise.md](20_plastid_canonicalise.md) —
  implement `bin/plastid_canonicalise.py` from the algorithm doc
  alone, without opening `combine_gfa.py`.

This task (18) ships C3 only, so the `plant_pt` branch emits the
dominant plastid contig (which may be LSC-only in the 3-edge case)
without splitting into `path1.fasta` / `path2.fasta` isoforms yet.
Task 20 tightens the plant_pt `target_min_bp` bound after `path1`
substitution lands.

**Prerequisite:** [task 16 — METAFLYE](completed/16_metaflye.md) is
real, so `BANDAGE_NG.out.assembly` (which forwards METAFLYE's
`assembly.fasta`, `assembly_graph.gfa`, `assembly_info.txt` + the
PNG from [task 17](17_bandage_ng.md)) carries real assembly artefacts
into this stage. Channel wiring already exists in
[main.nf:105](../main.nf#L105) —
`BIN_TARGET(BANDAGE_NG.out.assembly, ch_organelle_refs)` — and
requires no topology change.

**Exit criteria:**

- `-profile integration` produces, for every `status: "ok"` sample:
    - `results/<sample_id>/bin_target/target.fasta` — the selected
      dominant target contig (single-record FASTA for `animal_mt` and
      `plant_pt`; may be multi-record for `plant_mt`
      ([spec §3.3](../spec/03-organelles.md#33-specific-issues-and-decisions)
      multi-isoform).
    - `results/<sample_id>/bin_target/secondaries.tsv` — one row per
      non-selected contig with columns `contig_id`, `length_bp`,
      `coverage`, `ref_identity_pct`, `orf_ok`, `classification`
      (`off_target` | `secondary_target` | `low_confidence`).
    - `results/<sample_id>/bin_target/bin_metadata.json` — captures
      the decision rationale for the primary contig, including
      circularity for `animal_mt`.
- New integration assertion: `INT-ANIMAL-01/bin_target/target.fasta`
  contains **exactly one contig** in the range 14–20 kb (the pea
  aphid mitogenome
  [spec §3.1](../spec/03-organelles.md#31-what-differs-between-assembly-targets));
  `bin_metadata.json.circular == true`. `INT-PLANT-01-pt` yields a
  single contig in the range 100–170 kb (LSC+IR+SSC path, un-canonicalised
  is fine at this point).
- `BIN_TARGET` runs in a dedicated container that has Python 3.12 +
  `minimap2` + `biopython` on `PATH` (see §3).
- `bin/bin_target.py` has unit tests covering:
    - Selection when a single obvious contig dominates.
    - Selection under a coverage tie (identity + ORF break the tie).
    - Selection under an identity tie (coverage breaks the tie).
    - Non-selection when no contig clears the minimum thresholds
      (empty `target.fasta`; all rows in `secondaries.tsv`
      classified `low_confidence`).
    - Circularity check on a synthetic end-overlapping contig
      (positive case) and a linear contig (negative case).
- Both `-profile stub -stub-run` (fast CI) and `-profile integration`
  (nightly) go green end-to-end.

**Not in scope:**

- **C4 plastid canonicalisation.** `plastid_canonicalise.py` and its
  `path1.fasta` / `path2.fasta` outputs are deferred to a follow-up
  task; the module's commented `emit: isoforms` output at
  [modules/local/bin_target.nf:20-21](../modules/local/bin_target.nf#L20-L21)
  stays commented. Add the follow-up to
  [tasks/todo.md](todo.md) on completion.
- Plant `plant_mt` isoform classification. Multi-contig plant mt is
  ([spec §3.3](../spec/03-organelles.md#33-specific-issues-and-decisions))
  legitimate; ship all contigs passing the thresholds and let
  downstream stages consume the multi-record FASTA.
- Downstream conditional wiring on binning failure. If
  `target.fasta` is empty, `BLAST_VALIDATE` / `ANNOTATE` /
  `MINIPROT_EXTRACT` will produce empty outputs; `COLLATE` handles
  the `no_assembly` / `no_barcode` distinction
  ([spec principle 7](../CONSTITUTION.md)) as its own P4 task.
- Tunability of the coverage-spike / identity thresholds. Ship the
  defaults from §2 below; add a benchmarking follow-up if the
  integration fixtures push against them.

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pinned by SHA, never `latest`.
- No host tools; `minimap2`, `python3.12`, `biopython` all come from
  the container.
- `bin/bin_target.py` is auto-staged onto every process `PATH` by
  Nextflow's `bin/` convention — invoke it by bare name.
- Organelle refs enter via the existing `ch_organelle_refs` value
  channel — no bare `${params.organelle_refs}` in the script block.

---

## 1. Per-target genetic-code table — `nextflow.config`

Add a `genetic_code_tables` map alongside `assembly_size_hints` in
[`nextflow.config`](../nextflow.config). Sourced from
[spec §3.1](../spec/03-organelles.md#31-what-differs-between-assembly-targets):

```groovy
// NCBI genetic-code tables per assembly target — spec §3.1.
// animal_mt entry is a two-table clade trial per spec §3.3; C3 selects
// the table under which the primary contig yields the longest valid ORF.
genetic_code_tables = [
    animal_mt: [2, 5],   // vertebrate / invertebrate — trial both
    plant_pt:  [11],     // bacterial/plastid
    plant_mt:  [1],      // standard
]
```

The Python script consumes these as a comma-separated CLI flag (see
§2), not by reaching into `params` — keeps the script testable in
isolation and keeps Nextflow the single owner of parameter resolution.

**Thresholds** (defaults, hardcoded in `bin/bin_target.py` for P3;
promoted to `params.bin_target_thresholds` only if benchmarking shows
they need per-target tuning):

| Threshold | Default | Rationale |
|---|---|---|
| Minimum ref identity (%) | 80 | Coarse organelle match; RECRUIT already filtered off-kingdom. |
| Minimum ref-aligned fraction | 0.30 | ≥ 30 % of the contig aligns to the target reference — filters chimeras and long off-target survivors. |
| Coverage-spike factor | 2.0 | Selected contig's coverage ≥ 2× the median of all contigs (typical organelle-vs-nuclear ratio at skim depth). |
| Minimum ORF length (aa) | 100 | Below this, "ORF" is likely a random hit; irrelevant to true organelle protein-coding genes. |
| End-overlap window (`animal_mt`) | 300 bp | Flye typically leaves a ~overlap on closing a circle; 300 bp is generous enough to detect it without spurious matches. |

## 2. `bin/bin_target.py` (new — custom logic C3)

Python + `biopython` + `mappy` (Python bindings for minimap2). Uses
`biopython` for ORF translation across NCBI genetic-code tables and
`mappy` for in-process minimap2 alignment (no subprocess plumbing).

Pseudocode of the intended shape:

```python
#!/usr/bin/env python3
"""BIN_TARGET (C3) — spec §2 stage 10, §3.3.

Selects the dominant target contig from a METAFLYE assembly using the
intersection of:
  (a) coverage spike vs. per-contig median (from assembly_info.txt),
  (b) reference identity + aligned fraction to ${target}.mmi (minimap2),
  (c) ORF integrity under target-appropriate NCBI genetic-code table(s).

Emits target.fasta, secondaries.tsv, and bin_metadata.json describing
the decision. animal_mt: also runs the end-overlap circularity check.
"""

# --- CLI ---
# --assembly assembly.fasta
# --assembly-info assembly_info.txt
# --organelle-ref  <target>.mmi
# --sample-id      <id>
# --assembly-target {animal_mt, plant_pt, plant_mt}
# --genetic-codes  <csv of NCBI table ids>
# --out-target     target.fasta
# --out-secondaries secondaries.tsv
# --out-metadata    bin_metadata.json

# --- Algorithm ---
# 1. Load assembly_info.txt → per-contig (length, coverage).
# 2. For each contig: minimap2 (mappy) vs the organelle ref → (identity, aligned_frac).
# 3. For each contig: longest ORF (any frame, +/-) under each supplied
#    genetic-code table; record best-table + best-ORF-length.
# 4. Classify:
#      - target_candidate  : identity ≥ 80 AND aligned_frac ≥ 0.30
#                            AND coverage ≥ 2 × median coverage
#                            AND ORF ≥ 100 aa
#      - secondary_target  : passes (a) + (b) but not coverage spike
#      - off_target        : fails identity OR aligned_frac
#      - low_confidence    : passes some signals but not the intersection
# 5. Select primary:
#      - animal_mt / plant_pt: single best target_candidate by
#        (coverage_rank, identity_rank, orf_length_rank).
#      - plant_mt: all target_candidates (multi-contig legitimate).
# 6. animal_mt only: end-overlap circularity check on the primary
#    (mappy self-align first vs last 300 bp of the contig; ≥ 90 %
#    identity on ≥ 200 bp → circular=True).
# 7. Emit target.fasta, secondaries.tsv, bin_metadata.json.
# 8. Always exit 0 unless a tool crashes — an empty target.fasta is a
#    valid, informative outcome (COLLATE will surface it as no_assembly).
```

Rationale for the shape:

- **`mappy` not subprocess.** Avoids a `Popen` layer, keeps the code
  and unit tests direct. `mappy` is the official Python interface to
  minimap2 and is available on Bioconda.
- **`biopython` for ORF search.** `Bio.SeqUtils.six_frame_translations`
  or a direct `Seq.translate(table=N)` loop is short, tested, and
  correct across NCBI tables — no need to hand-roll.
- **Selection function is pure.** `select_primary(rows) -> (primary, secondaries)`
  takes the enriched per-contig rows and returns the split; unit
  tests target this function directly with synthetic rows, no I/O.
- **Always exit 0 on decisions.** An empty `target.fasta` is data, not
  a Nextflow error — parallel to
  [COVERAGE_GATE](completed/15_coverage_gate.md)'s soft-fail contract.

## 3. Container build

Use `mulled-build` (as in
[task 15 §3](completed/15_coverage_gate.md#3-container-build)) to produce
a conda-based mulled image with Python 3.12 + `minimap2` + `biopython`
+ `mappy`:

```bash
mulled-build build 'python=3.12,minimap2=2.31,biopython=1.85,mappy=2.31'
```

`mulled-build` prints the generated `quay.io/biocontainers/mulled-v2-…`
tag. Retag with the tool-version convention used for other C-component
containers and push:

```bash
docker tag \
    quay.io/biocontainers/mulled-v2-<hash>:<hash>-0 \
    neoformit/daff-wf5-bin-target:python3.12_minimap2-2.31_biopython-1.85
docker push neoformit/daff-wf5-bin-target:python3.12_minimap2-2.31_biopython-1.85
```

Pin by digest in [`conf/containers.config`](../conf/containers.config),
removing the `TODO P3` comment at
[conf/containers.config:23-26](../conf/containers.config#L23-L26):

```groovy
withName: 'BIN_TARGET' {
    container = 'neoformit/daff-wf5-bin-target:python3.12_minimap2-2.31_biopython-1.85@sha256:<pinned>'
}
```

**Note on C3/C4 shared image.** [Spec §2.2](../spec/02-stages.md#22-custom-logic-components)
calls for C3 and C4 to share the image `wf5/bin-target:<tag>` because
C4 is invoked as a helper by C3. This task builds the shared image
even though C4 is not yet wired — the image only needs Python +
`biopython` + `minimap2`/`mappy`, and C4 (GFA parsing + sequence
concatenation) needs no additional deps. When [task N — C4 plastid
canonicalisation](todo.md) lands, it will drop `plastid_canonicalise.py`
into `bin/` and consume the same container without a rebuild.

## 4. `modules/local/bin_target.nf` (replace stub)

```groovy
// Stage 10 — per-contig binning: coverage spike ∩ ref identity ∩ ORF integrity.
// C3 custom logic — see spec §2.2, §3.3.
// Plant-cp canonicalisation (C4) and its path1/path2 isoform outputs
// deferred to a follow-up task; the commented emit block below stays
// commented until then.
// Animal-mt: end-overlap circularity check recorded in bin_metadata.json.

process BIN_TARGET {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/bin_target",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(assembly), path(gfa), path(info), path(graph_png)
    path organelle_refs

    output:
    tuple val(meta), path("target.fasta"), path("secondaries.tsv"), emit: binned
    path("bin_metadata.json"), emit: metadata
    // TODO(task-20 C4): uncomment when plastid_canonicalise.py lands.
    // path "plastid_isoforms/", optional: true, emit: isoforms

    stub:
    """
    touch target.fasta secondaries.tsv bin_metadata.json
    """

    script:
    def codes = params.genetic_code_tables[meta.assembly_target].join(',')
    """
    bin_target.py \\
        --assembly ${assembly} \\
        --assembly-info ${info} \\
        --organelle-ref ${organelle_refs}/${meta.assembly_target}.mmi \\
        --sample-id ${meta.sample_id} \\
        --assembly-target ${meta.assembly_target} \\
        --genetic-codes ${codes} \\
        --out-target target.fasta \\
        --out-secondaries secondaries.tsv \\
        --out-metadata bin_metadata.json
    """
}
```

Notes:

- **Output tuple shape unchanged for `binned`.** Downstream consumers
  at [main.nf:108, 124](../main.nf#L108) already expect
  `(meta, target_fasta, secondaries)`. The new `metadata` emit is a
  separate channel to be picked up by COLLATE when P4 lands; the
  current [main.nf:124-142](../main.nf#L124-L142) join is not
  touched here.
- **Reference path pattern.** `ch_organelle_refs` is the bundle root
  directory; the script indexes into it as
  `${organelle_refs}/${meta.assembly_target}.mmi`. This mirrors the
  pattern used by [`RECRUIT`](../modules/local/recruit.nf) and
  keeps the reference layout single-source-of-truth
  ([spec §4.4](../spec/04-reference-data.md#44-consolidated-build-script)).
- **`process_medium` label.** Alignment + ORF search on organelle-
  scale assemblies is modest; keep at medium unless benchmarking on
  a very fragmented `plant_mt` fixture pushes it up.
- **No `set -euo pipefail`.** The Python script owns exit-code
  semantics and always returns 0 on a valid decision.

## 5. Unit tests — `scripts/tests/test_bin_target.py`

Independent of Nextflow. Fixtures constructed inline in Python — no
checked-in binary test data (drift is a maintenance burden per
[constitution rule 19](../CONSTITUTION.md)).

Cases:

| # | Scenario | Expected outcome |
|---|---|---|
| 1 | 1 high-cov high-identity ORF-passing contig + 3 low-cov off-target survivors | Primary = that contig; 3 secondaries all `off_target` |
| 2 | 2 high-cov contigs, only one hits ORF | Primary = ORF-passing one; other `low_confidence` |
| 3 | Coverage tie (2 contigs at identical cov, both hit ref) | Tie-breaker rules: identity > ORF length; assert stable primary |
| 4 | No contig clears identity | Empty `target.fasta`; all rows `off_target` or `low_confidence`; exit 0 |
| 5 | animal_mt synthetic circular contig (300 bp end-overlap) | `bin_metadata.json.circular == true` |
| 6 | animal_mt synthetic linear contig | `bin_metadata.json.circular == false` |
| 7 | plant_mt multi-contig valid input (3 contigs all pass) | `target.fasta` contains all 3 records |
| 8 | animal_mt genetic-code trial: contig has valid ORF only under table 5 | `bin_metadata.json.selected_genetic_code == 5` |

Fixture construction: build synthetic FASTA + `assembly_info.txt`
inline; use a tiny handcrafted `.mmi` built from a fragment of the
GetOrganelle `animal_mt` panel (small enough to check in as a
`scripts/tests/fixtures/mini_animal_mt.mmi`, < 100 kB) so `mappy`
alignment is deterministic.

**Coverage target:** 100 % branch coverage of the selection function
per [constitution rule 14](../CONSTITUTION.md) — C3 is a custom-logic
component.

## 6. Integration-test wiring

Add a BIN_TARGET assertion block after the METAFLYE block in
[`tests/integration/assertions.sh`](../tests/integration/assertions.sh):

```bash
# BIN_TARGET is real (task 18):
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    tgt="$OUTDIR/$sample/bin_target/target.fasta"
    meta="$OUTDIR/$sample/bin_target/bin_metadata.json"
    if [[ ! -s "$tgt" ]]; then
        echo "FAIL: $sample target.fasta missing or empty"
        FAILED=1
        continue
    fi
    n_contigs=$(grep -c '^>' "$tgt")
    tgt_bp=$(grep -v '^>' "$tgt" | tr -d '\n' | wc -c)
    target="${SAMPLE_TARGET[$sample]}"
    bounds="tests/integration/expected/${target}/bin_bounds.json"
    min_bp=$(jq .target_min_bp "$bounds")
    max_bp=$(jq .target_max_bp "$bounds")
    max_contigs=$(jq .target_max_contigs "$bounds")
    if (( tgt_bp < min_bp || tgt_bp > max_bp )); then
        echo "FAIL: $sample target ${tgt_bp} bp outside [${min_bp}, ${max_bp}]"
        FAILED=1
    elif (( n_contigs > max_contigs )); then
        echo "FAIL: $sample $n_contigs contigs > max ${max_contigs}"
        FAILED=1
    else
        echo "OK:   $sample target ${tgt_bp} bp / ${n_contigs} contigs"
    fi
    # animal_mt: additionally assert circularity.
    if [[ "$target" == "animal_mt" ]]; then
        circular=$(jq -r .circular "$meta")
        if [[ "$circular" != "true" ]]; then
            echo "FAIL: $sample animal_mt target not detected circular"
            FAILED=1
        fi
    fi
done
```

**New fixtures required** — one per target under
[`tests/integration/expected/`](../tests/integration/expected/):

`tests/integration/expected/animal_mt/bin_bounds.json`:
```json
{
    "target_min_bp": 14000,
    "target_max_bp": 20000,
    "target_max_contigs": 1,
    "circular_expected": true
}
```

`tests/integration/expected/plant_pt/bin_bounds.json`:
```json
{
    "target_min_bp": 70000,
    "target_max_bp": 170000,
    "target_max_contigs": 3
}
```

> **Note on the `plant_pt` lower bound.** Pre-canonicalisation,
> `target.fasta` on `plant_pt` may be the LSC edge alone (~85 kb)
> when Flye emits the plastid as 3 edges rather than a single
> resolved circle. The loose 70 kb floor accepts either the LSC-only
> outcome or a resolved circle. Task 20 (C4) tightens this to
> ~140 kb once `path1` substitution lands.

`tests/integration/expected/plant_mt/bin_bounds.json`:
```json
{
    "target_min_bp": 200000,
    "target_max_bp": 700000,
    "target_max_contigs": 20
}
```

Bounds are **coarser than the METAFLYE `assembly_bounds.json`** on
purpose: METAFLYE bounds are total-assembly (includes off-target
survivors under `--meta`); BIN_TARGET bounds are the target contig(s)
only, so the lower bound is tighter and the upper bound is closer to
the true organelle size. See
[task 16 §10 outcomes](completed/16_metaflye.md#10-outcomes) for the
observed animal_mt case (15 contigs / 349 kb total → 1 circular
~17 kb contig after BIN_TARGET).

Also update the progressive-uncomment header at the top of
`assertions.sh` to tick off BIN_TARGET.

## 7. Fast CI

No change. Stub script block still `touch`es all three outputs
(`target.fasta`, `secondaries.tsv`, `bin_metadata.json`);
container-coverage lint stays satisfied.

## 8. Follow-up tasks

Two sibling tasks under a clean-room-firewall design (upstream ptGAUL
`combine_gfa.py` is unlicensed):

1. [19_plastid_algorithm_spec.md](19_plastid_algorithm_spec.md) —
   read the ptGAUL source, produce `spec/plastid-canonicalisation.md`.
2. [20_plastid_canonicalise.md](20_plastid_canonicalise.md) —
   implement `bin/plastid_canonicalise.py` from that spec alone,
   wire the `emit: isoforms` output of `BIN_TARGET`, and apply the
   `path1`-for-`target.fasta` substitution that tightens the
   `plant_pt` bounds.

Task 20 reuses the `neoformit/daff-wf5-bin-target` image built in
this task; no container work in either 19 or 20. No `todo.md` entry
needed.

## 9. Verification

1. `pytest scripts/tests/test_bin_target.py` — all 8 cases pass;
   coverage report shows 100 % branch on the selection function.
2. `flake8 bin/bin_target.py scripts/tests/test_bin_target.py
   --max-line-length=100` — clean.
3. `nextflow run . -profile stub -stub-run` — 17/17 processes hit
   their stub, no regressions.
4. `nextflow run . -profile integration` locally:
   - `INT-ANIMAL-01/bin_target/target.fasta` → 1 contig, ~17 kb,
     `bin_metadata.json.circular == true`,
     `selected_genetic_code == 5` (invertebrate).
   - `INT-PLANT-01-pt/bin_target/target.fasta` → 1 contig, ~150 kb,
     `selected_genetic_code == 11`.
   - `INT-PLANT-01-mt/bin_target/target.fasta` → 1–5 contigs, total
     ~250–700 kb, `selected_genetic_code == 1`.
5. `bash tests/integration/assertions.sh` — new BIN_TARGET block
   reports OK for all three samples; pre-existing blocks still green.
6. CI green on `Tests` (fast) + nightly `Integration` on the PR.

## 10. Deliverables checklist

- [x] `bin/bin_target.py` — real implementation, `mappy` + `biopython`,
      always-exit-0 on decisions.
- [x] `scripts/tests/test_bin_target.py` — 10 cases (8 spec cases + 2
      additional edge cases); 100 % branch coverage on the selection function.
      Tests use in-memory synthetic sequences — no binary fixtures required.
- [x] Container built via `mulled-build`, retagged as
      `neoformit/daff-wf5-bin-target:python3.12_minimap2-2.31_biopython-1.85`,
      pushed to Docker Hub.
- [x] [`conf/containers.config`](../conf/containers.config) — SHA-pinned
      container for `BIN_TARGET`, `TODO P3` comment removed.
- [x] [`nextflow.config`](../nextflow.config) —
      `params.genetic_code_tables` map added.
- [x] [`modules/local/bin_target.nf`](../modules/local/bin_target.nf) —
      real `script:` block invoking `bin_target.py`; stub retained;
      new `metadata` emit added; C4 `isoforms` emit stays commented.
- [x] [`tests/integration/expected/*/bin_bounds.json`](../tests/integration/expected/) —
      three new fixtures.
- [x] [`tests/integration/assertions.sh`](../tests/integration/assertions.sh) —
      BIN_TARGET assertion block + circularity check for animal_mt;
      progressive-uncomment header updated.
- [ ] Fast CI (`Tests`) + `Integration` workflows both green on the PR.

## 11. Notes / non-issues

- **No channel-topology change to the `binned` emit.** The existing
  `main.nf:124` join keys off `(meta, target_fasta, secondaries)`
  and continues to work. The new `metadata` emit is picked up by
  `COLLATE` in P4 — not this task.
- **Empty `target.fasta` is a valid outcome.** Downstream stages will
  produce empty outputs; `COLLATE`'s P4 task will introduce the
  `no_assembly` bundle path
  ([spec principle 7](../CONSTITUTION.md)).
- **`plant_mt` multi-record output.** Legitimate per
  [spec §3.3](../spec/03-organelles.md#33-specific-issues-and-decisions);
  the `max_contigs: 20` bound in the fixture is generous — tighten
  only if the SRR11315861 fixture consistently returns fewer.
- **Circularity method choice.** End-overlap detection is the
  cheapest signal; a full DFS on the GFA is a P3+ refinement if the
  end-overlap check misfires on tangled assemblies. Not this task.
- **Genetic-code table selection is recorded, not enforced.** C5
  (`validate_barcodes.py`, task TBD) reruns the clade trial per-locus
  in `MINIPROT_EXTRACT`; C3's choice is a diagnostic for reporting
  ([spec §3.3](../spec/03-organelles.md#33-specific-issues-and-decisions))
  and can safely differ from C5's per-locus decision.

## 12. Outcomes

- Container: `neoformit/daff-wf5-bin-target:python3.12_minimap2-2.31_biopython-1.85@sha256:532a372c92e6e3152be9f1da4c2f4201e44d76721900fd8faa03cf2bebce1698`
  Built from `mulled-build 'python=3.12,minimap2=2.31,biopython=1.85,mappy=2.31'`.
- Selection function branch coverage: 100% on `select_primary`; 10 tests all passing.
- Unit tests use in-memory synthetic sequences only (no `.mmi` fixture checked in).
  Case 8 tests table 5 vs table 2 via the AGA codon difference (stop in t2, Ser in t5),
  exercised through the RC reading frame — the RC of TCT is AGA.
- Integration run pending (CI).
- Any threshold tuning applied: none.
