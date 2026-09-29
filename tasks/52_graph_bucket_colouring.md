# Task 52 — Move `BANDAGE_NG` after `BIN_TARGET`; colour the assembly graph by target/secondary/off-target

**Phase:** P4a (spec 06-phases.md).
**Goal:** Reorder stage 9 (`BANDAGE_NG` + C12, the node-labelled assembly
graph SVG from task 47_bandage_interactive_svg.md) to run after stage 10
(`BIN_TARGET`), and use the binning classification that becomes
available there to colour each graph node by which bucket(s) — target,
secondary, off-target — traverse it, instead of the flat neutral grey
task 47 shipped.

**Depends on:** task 47_bandage_interactive_svg.md (ships the graph SVG,
the sentinel round-trip, and `data-contigs` per node — this task adds
colour on top, it does not change the round-trip mechanism).

**Related tasks:** task 53_graph_bootstrap_tooltips.md (Bootstrap-styled
hover tooltips for the same graph) is independent of this one — either
can land first — but both touch `bin/annotate_graph_svg.py` and the
graph's rendered `<g>` markup, so whichever lands second should re-read
the other's diff before starting.

---

## 1. Overview

Task 47 gave the assembly-graph diagnostic node identity (`data-node`,
`data-contigs`, a hover tooltip) but deliberately left every node the
same neutral grey. That was not an oversight — at the point stage 9
(`BANDAGE_NG`) runs in the pipeline, stage 10 (`BIN_TARGET`, the
component that decides which contig is the target organelle, which
are secondary candidates, and which are off-target) has not run yet.
The graph is built before the pipeline knows which contigs matter, so
task 47 §7 explicitly ruled out inventing a colour it couldn't yet
justify (CONSTITUTION principle 7) and flagged this as a follow-up:
*"probably colour by the set of buckets traversed, with a distinct
colour for mixed [...] it would have to move the annotation step after
stage 10."*

This task does exactly that. Once `BIN_TARGET` has run,
`bin_metadata.json` records, for every contig, which bucket it fell
into — the same three buckets (`target` / `secondary` / `off-target`)
the Assembly tab's per-contig coverage chart already colours
(`bin/report/report.py`'s `_contig_bucket()` / `_BUCKET_COLOURS`,
green/orange/grey). A GFA segment (`edge_N`) already carries, via
task 47's `data-contigs`, the set of contigs whose Flye `graph_path`
traverses it. Joining that set against `bin_metadata.json`'s
per-contig classification gives every segment a bucket — usually one,
sometimes more than one, since segments and contigs are a many-to-many
relationship (spec/02-stages.md's C12 row, task 47 §2 item 6): a
plastid inverted-repeat edge is a real example of a segment that can
be shared between the target contig and an off-target one.

**Reordering is the enabling change, not a big one on its own.**
`BANDAGE_NG`/C12 only ever needed `assembly_graph.gfa` and
`assembly_info.txt` (both from `METAFLYE`, stage 7) to do its job;
today's dependency on running before `BIN_TARGET` is purely
positional, not a real data dependency — `bin_target.py` never reads
the graph SVG or the sentinel CSV it is currently (needlessly) handed
through `modules/local/bin_target.nf`'s `graph_svg` input. Swapping the
two stages' order removes that unused pass-through parameter as a
side effect.

### Why "after `BIN_TARGET`", not "at the end of the run"

The prompt that raised this idea suggested moving the graph stage
"near the end of the run." The binning classification is finalised at
stage 10 and nothing downstream of it (BLAST validation, CDS
annotation, barcode extraction) refines *which bucket a contig is in*
— those stages consume the binning decision, they don't revise it. So
running the graph any later than immediately after `BIN_TARGET` buys
nothing for the colouring this task adds, while giving something up:
today `BANDAGE_NG` runs in parallel with the independent stages 11–14
(`BLAST_VALIDATE`, `MINIPROT_CDS`, barcode extraction, annotation), and
its diagnostic is available for a human to look at while those still
run. Pushing it to the very end would serialise it behind work it has
no dependency on and delay the one output most useful for spotting a
messy assembly early. Recommendation: swap stages 9 and 10 only.

---

## 2. What the current code already gives this task

Verified by reading, not re-derived here:

- `bin/report/report.py`'s `_contig_bucket(classification)` already
  maps `bin_metadata.json`'s per-contig `classification` field to one
  of three buckets: `target_candidate` → `target`, `secondary_target`
  → `secondary`, anything else (`off_target`, `sibling_organelle`) →
  `off-target`. `_BUCKET_COLOURS` already pins the palette:
  `target` `#2ca02c` (green), `secondary` `#ff7f0e` (orange),
  `off-target` `#7f7f7f` (grey) — the same colours the per-contig
  coverage chart on the Assembly tab uses today. Reuse these values
  for the graph so the two pictures agree (rule 19); do not invent a
  second palette.
- `bin/annotate_graph_svg.py`'s `parse_graph_path_column()` (task 47)
  already returns `{segment_name: {contig_id, ...}}` — exactly the
  join key this task needs against `bin_metadata.json`'s per-contig
  classifications.
- `modules/local/bin_target.nf`'s `output: ... emit: metadata` already
  publishes `bin_metadata.json` per sample, joinable by `meta` — no
  new Nextflow output needed, just a new input wired to
  `ANNOTATE_GRAPH_SVG`.

---

## 3. Design

### 3.1 Stage renumbering

Swap stage 9 and stage 10 throughout the spec and code comments:
`BIN_TARGET` becomes stage 9, `BANDAGE_NG` (+ C12) becomes stage 10.
No other stage number changes. Update:

- spec/01-pipeline-flow.md's flow diagram (§1).
- spec/02-stages.md's stage table (§2) and the C12 row in §2.2 (also
  update its container-choice comment, which currently reconciles the
  three-process split against task 47 §4.2 — the process split itself
  is unchanged by this task, only its position in the DAG).
- spec/06a-reports.md's Assembly tab row (§3) — the "node colour
  carries no binning meaning" sentence is no longer true and must be
  rewritten to describe what colour *does* mean now (§3.3 below).
- `modules/local/bandage_ng.nf`, `modules/local/allocate_graph_sentinels.nf`,
  `modules/local/annotate_graph_svg.nf`, `modules/local/bin_target.nf`:
  their "Stage 9"/"Stage 10" header comments.
- `bin/README.md`'s component table (C3's stage 10 → 9, C12's stage
  9 → 10).
- `tests/integration/assertions.sh`'s progressive-uncomment-plan
  comment block at the top (§ "BANDAGE_NG (real renderer)" /
  "BIN_TARGET (real C3)" rows).

### 3.2 `main.nf` reordering

```
METAFLYE
  -> BIN_TARGET(METAFLYE.out.assembly, ch_organelle_refs)
       # consumes (meta, assembly, gfa, info) directly again —
       # the graph_svg pass-through parameter this task removes
  -> ALLOCATE_GRAPH_SENTINELS(METAFLYE.out.assembly)
       # unchanged: still only needs the GFA
  -> BANDAGE_NG(ALLOCATE_GRAPH_SENTINELS.out.sentinels)
       # unchanged
  -> ANNOTATE_GRAPH_SVG(
       BANDAGE_NG.out.rendered.join(BIN_TARGET.out.metadata, by: 0))
       # new: bin_metadata.json joined in as an extra input
  -> BLAST_VALIDATE(BIN_TARGET.out.binned)
       # unchanged, no longer waits on the graph chain
```

`ch_graph_svg` (built from `ANNOTATE_GRAPH_SVG.out.assembly` today)
keeps the same shape for `COLLATE`; only where it originates in the
DAG changes.

### 3.3 `modules/local/annotate_graph_svg.nf`

Add `path(bin_metadata_json)` to the input tuple (from
`BIN_TARGET.out.metadata`, joined by `meta`) and pass it to
`annotate_graph_svg.py annotate` as a new `--bin-metadata` argument,
mirroring how `bin/organelle_map.py` already takes the same file
(`load_bin_metadata()` — task 44) as an optional path that degrades to
"nothing known" when absent/malformed, not an error.

### 3.4 `modules/local/bin_target.nf`

Drop the now-unused `path(graph_svg)` from the input tuple — it was
never read by `bin_target.py`, only carried through positionally.

### 3.5 `bin/annotate_graph_svg.py` — bucket colouring

Factor the classification → bucket → colour logic that
`bin/report/report.py` already owns (`_contig_bucket`,
`_BUCKET_COLOURS`) out into a small shared module — e.g.
`bin/contig_bucket.py`, exporting `contig_bucket(classification) ->
str` and `BUCKET_COLOURS: dict` — imported by both `report.py` and
`annotate_graph_svg.py` (rule 19: one definition, not two copies that
can drift). `annotation_gff.py`, already shared between C10/C11
(task 40), is the precedent for this kind of small shared parsing/logic
module living in `bin/` without its own Nextflow process.

New function in `annotate_graph_svg.py`, used by `annotate_svg()` in
place of the current always-neutral fill:

```
def segment_bucket(contig_ids: set, contig_classifications: dict) -> str | None:
    """contig_classifications: {contig_id: bin_metadata classification}.
    Returns one of BUCKET_COLOURS' keys, or 'mixed' if the segment's
    traversing contigs span more than one bucket, or None if no
    traversing contig's classification is known (falls back to the
    existing neutral display colour)."""
    buckets = {
        contig_bucket(contig_classifications[c])
        for c in contig_ids if c in contig_classifications
    }
    if not buckets:
        return None
    if len(buckets) > 1:
        return 'mixed'
    return next(iter(buckets))
```

A `mixed` segment gets a fourth, distinct colour not used elsewhere in
the report (e.g. a blue/purple — pick something that reads clearly
against the existing green/orange/grey trio and the white background;
exact hex is an implementation choice, just keep it out of
`BUCKET_COLOURS`' own palette so it can't be confused with a real
bucket). A segment with no known bucket (traversing contig(s) present
in `data-contigs` but absent from `bin_metadata.json`, or no
traversing contig at all) keeps today's neutral
`DISPLAY_COLOUR` (`#c8c8c8`) — same graceful-degradation shape task 47
already uses for a missing `assembly_info.txt`.

`bin_metadata.json` absent or malformed (mirrors task 44's
`load_bin_metadata()`): every segment falls back to the neutral
colour, same as today — this must never turn into a hard failure
(CONSTITUTION principle 8, rule 8 — C12 already always exits 0).

### 3.6 Report template + info badge

`scripts/report/templates/components/assembly.html`'s graph info-badge
text currently says colour carries no binning meaning — rewrite it to
say what the four colours (green/orange/grey/mixed) mean, matching the
coverage-chart legend's wording so a reader doesn't have to reconcile
two descriptions of the same palette. Add a small colour-key legend
near the graph if the existing prose isn't enough on its own (the
coverage chart already has its own Plotly-generated legend for
comparison — check whether reusing that pattern or a plain HTML swatch
list reads better here).

---

## 4. Work items

1. Swap the stage-9/stage-10 labels and reorder the Nextflow DAG per
   §3.1–3.2.
2. `bin/contig_bucket.py` (new, shared) — extract `contig_bucket()` /
   `BUCKET_COLOURS` from `bin/report/report.py`; update `report.py` to
   import from it instead of defining its own copy.
3. `bin/annotate_graph_svg.py` — add `--bin-metadata` to the
   `annotate` subcommand, `segment_bucket()`, and wire it into
   `annotate_svg()`'s fill-colour choice (replacing the current
   always-`DISPLAY_COLOUR` behaviour for segments with a known
   bucket).
4. `modules/local/annotate_graph_svg.nf` / `modules/local/bin_target.nf`
   / `main.nf` per §3.2–3.4.
5. `scripts/report/templates/components/assembly.html` info-badge +
   legend per §3.6.
6. Spec + comment reconciliation per §3.1.

---

## 5. Tests

### 5.1 Unit — `scripts/tests/test_contig_bucket.py` (new)

100% branch coverage (rule 14): every `classification` value the real
pipeline emits (`target_candidate`, `secondary_target`, `off_target`,
`sibling_organelle`) plus an unrecognised/`None` value.

### 5.2 Unit — `scripts/tests/test_annotate_graph_svg.py` (extend)

- `segment_bucket()`: single-bucket cases (target/secondary/
  off-target); mixed (≥2 distinct buckets, including the 3-way case);
  a contig id present in `data-contigs` but absent from
  `bin_metadata.json`; empty `contig_ids`.
- `annotate_svg()` / `run_annotate()`: a segment coloured by its
  bucket; a `mixed` segment gets the distinct fourth colour, not a
  bucket colour or the neutral one; missing/malformed
  `--bin-metadata` falls back to the neutral colour for every segment,
  same as today, and does not raise.

### 5.3 Unit — `scripts/tests/test_report.py` (extend)

Update `_contig_bucket`/`_BUCKET_COLOURS` references to import from
`bin/contig_bucket.py`; existing bucket-colour assertions for the
coverage chart must keep passing unchanged (this task must not change
what the chart shows, only add the same information to the graph).

### 5.4 Integration — `tests/integration/assertions.sh`

Extend the existing `graph.svg` structural assertion (task 47 §5.3):
for `INT-ANIMAL-01`, assert the target contig's colour
(`#2ca02c`) appears at least once among the graph's `data-node`
elements — `bin_metadata.json`'s `contigs_selected` names the target
contig directly, so the assertion can look up which segment(s) that
contig's `graph_path` touches without hardcoding an edge number.

### 5.5 Fixtures

No new fixture data needed — the existing `INT-ANIMAL-01` /
`INT-PLANT-01-pt` / `INT-PLANT-01-mt` fixtures already produce
`bin_metadata.json` with real classifications. A fresh
`-profile integration` run is required (stage reordering changes the
DAG; task 47's own re-run obligation pattern applies here for the same
reason) — confirm `tests/integration/assertions.sh` stays fully green
and hand-check that at least one rendered `report.html` shows a
non-grey graph node.

---

## 6. Explicit non-goals

- Changing what `BIN_TARGET`'s classification means, or its
  thresholds (`params.bin_target_thresholds`) — this task only
  *displays* the existing decision on a second picture.
- A legend/interaction more elaborate than a static colour key (no
  click-to-filter-by-bucket, no toggling buckets on/off) — out of
  scope, same boundary task 47 §9 already drew for the graph.
- Colouring the organelle map (task 44_organelle_map.md) by bucket —
  it already colours by feature kind (CDS/tRNA/rRNA), a different and
  unrelated axis; not touched here.

---

## 7. Acceptance criteria

1. `BANDAGE_NG` (+ C12) runs after `BIN_TARGET` in the Nextflow DAG;
   stage numbers are swapped consistently across spec, module
   comments and `bin/README.md`.
2. A graph node traversed only by the target contig renders
   `#2ca02c`; only-secondary renders `#ff7f0e`; only-off-target (or
   unknown) renders the existing neutral grey; a node traversed by
   contigs in more than one bucket renders a fourth, distinct colour.
3. `bin_metadata.json` missing or malformed does not fail the sample
   (exit 0) and every node falls back to today's neutral colour.
4. `bin/contig_bucket.py` is the single definition of the
   classification → bucket → colour mapping; `report.py`'s coverage
   chart and the graph agree on the same three bucket colours.
5. New/changed `bin/*.py` reach 100% branch coverage under
   `scripts/pytest.sh`; `flake8` clean.
6. `nextflow run . -profile stub -stub-run` green with the reordered
   DAG.
7. A fresh `-profile integration` run completes;
   `tests/integration/assertions.sh` is fully green, including the
   extended bucket-colour assertion.
8. Info badge / legend text on the Assembly tab accurately describes
   what graph node colour now means; no page claims colour is
   meaningless anywhere it still appears (spec + template).
