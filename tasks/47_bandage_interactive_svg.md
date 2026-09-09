# Task 47 — `BANDAGE_NG`: node-labelled, hover-interactive assembly graph SVG

**Phase:** P4a (spec 06-phases.md).
**Goal:** Replace stage 9's flat graph PNG with an SVG whose nodes carry
their GFA segment names, so the Assembly tab's graph can name what the
reader is hovering over instead of being an un-navigable picture.

**Depends on:** task 43b_report_stage_tabs.md (shipped the Assembly tab
that displays the graph) and task 46_report_enhancement.md (moved the
graph into the two-column "Genome annotation" section it will occupy).

**Related tasks:** task 44_organelle_map.md owns the *other* inline SVG
in the same section; the two are independent but land in adjacent
markup, so whichever goes second should re-check the row layout.

---

## 1. Overview

Stage 9 runs BandageNG over Flye's `assembly_graph.gfa` and saves a
picture of the assembly graph. The graph is the honest picture of what
the assembler actually had to work with: each node is a stretch of
sequence, each connection is a place where the assembler could not
decide how the pieces join. A clean organelle assembly looks like one
simple circle; a messy or contaminated one looks like a tangle. This is
why it is a **diagnostic** — an operator looks at it to understand *why*
an assembly came out the way it did.

The problem is that today it ships as a **PNG**: a flat raster with no
information in it beyond the pixels. A reader can see there are four
blobs and a tangle, but cannot tell *which* node is which, so the
picture cannot be connected to any of the numbers in the per-contig
table beside it. To answer "which of these is my target contig?" the
operator has to leave the report, open BandageNG locally, and load the
GFA by hand.

This task makes the graph self-describing: hovering a node tells you its
segment name (`edge_3`) and which contig(s) run through it.

The obvious implementation — "ask BandageNG for SVG instead of PNG" —
**does not work on its own**, and §2 records why in detail, because it
is the single fact this task's design turns on. BandageNG *can* write
SVG, but the SVG it writes has no node identity in it whatsoever. The
identity has to be smuggled through the one channel BandageNG does
expose: node colour. §3 specifies that round-trip.

Spec 02-stages.md stage 9, spec 02-stages.md §2.2 (component C12), spec
01-pipeline-flow.md, spec 06a-reports.md and spec 00-overview.md have
**already been updated** to describe this design; this task implements
what they now say.

---

## 2. What was verified before this brief was written

Run against the pinned container
`quay.io/biocontainers/bandage_ng:2026.6.1--hca0ed12_0` and the real
`tests/integration/output/INT-PLANT-01-pt/assembly/assembly_graph.gfa`
fixture (3 segments: `edge_1`, `edge_2`, `edge_3`). Do not re-litigate
these; they are settled.

1. **`BandageNG image` accepts an `.svg` output path.** The format
   change itself is available — the help text names `.jpg`, `.png` and
   `.svg`.

2. **The SVG carries no node identity.** It is a Qt painter dump —
   `<desc>Generated with Qt</desc>` — containing 10 anonymous `<path>`
   elements for a 3-segment graph. **Zero `id` attributes, zero
   `<text>` elements.** Paths do not even correspond 1:1 to nodes, so
   there is nothing to key a tooltip off and no label text to read.

3. **There is no flag that adds labels.** The `image` subcommand's
   entire option set is `--height`, `--width`, `--color`. There is no
   `--names` / `--labels`.

4. **`--color` is a usable identity channel.** Given a CSV assigning
   each segment a unique colour, each node renders as **exactly one
   path with that exact fill**. The full colour inventory of the output
   was:

   ```
   1 fill="#ff0001"    ← edge_1
   1 fill="#ff0002"    ← edge_2
   1 fill="#ff0003"    ← edge_3
   1 fill="#ffffff"    ← background
   27 stroke="#000000" ← edge lines / outlines
   ```

   Exact string match, one fill per node, no collisions with the
   background or the link strokes. SVG (unlike PNG) preserves the
   colour as a literal attribute value, so there is no anti-aliasing
   or quantisation risk.

5. **The colour CSV requires a header row.** A headerless file is
   rejected with `colors.csv didn't contain color`. The header used was
   `name,color`.

6. **Segments and contigs are different namespaces, joined
   many-to-many.** Flye's GFA segments are `edge_N`; its contigs are
   `contig_N`; `assembly_info.txt`'s `graph_path` column maps between
   them, e.g. `contig_1 → -2,1,2` (traverses `edge_2` reversed, then
   `edge_1`, then `edge_2`). One contig spans several edges, and one
   edge can be traversed by several contigs. **This is why this task
   does not colour nodes by target/secondary/off-target bucket** — a
   repeat edge shared between a target and an off-target contig has no
   single correct colour. See §7.

---

## 3. The design — sentinel-colour round-trip

Three steps inside the one `BANDAGE_NG` process:

1. **Allocate sentinels.** Read the segment names from
   `assembly_graph.gfa` (`S` records, column 2). Assign each a unique
   `#RRGGBB` and write the `name,color` CSV.
2. **Render.** `BandageNG image ... --color sentinels.csv`, output
   `<sample_id>.graph.svg`.
3. **Annotate.** Post-process the SVG: for each `fill="#rrggbb"` that
   matches a sentinel, replace the fill with the neutral display colour
   and attach the segment's identity.

Steps 1 and 3 are component **C12** (`bin/annotate_graph_svg.py`).

**Sentinel allocation must be deterministic** — derived from segment
order, not random — so that re-running the pipeline on the same input
produces a byte-identical SVG. A non-deterministic colour would make
the output churn under `-resume` and defeat reproducibility
(CONSTITUTION rule 17).

Sentinel colours are a **machine-readable tag only**. They are consumed
and overwritten by step 3 and must never survive into the published
SVG, and nothing downstream may read meaning into them. Pick a range
that cannot be confused with Bandage's own background (`#ffffff`) or
link strokes (`#000000`).

### 3.1 What C12 attaches to each node path

- `data-node="edge_3"` — the segment name.
- `data-contigs="contig_1,contig_4"` — the contig(s) whose `graph_path`
  traverses this segment, comma-separated; empty attribute if none do.
- A `<title>` child element — this is what makes the basic hover label
  work **with no JavaScript at all**, since browsers render `<title>`
  inside SVG as a native tooltip. Prefer this over a JS tooltip: the
  report is a self-contained file that must degrade gracefully.

Example of the intended transformation (illustrative, not literal
output):

```
before:  <path fill="#ff0003" d="..."/>

after:   <g data-node="edge_3" data-contigs="contig_1">
           <title>edge_3 — contig_1</title>
           <path fill="#c8c8c8" d="..."/>
         </g>
```

### 3.2 Failure behaviour

C12 **always exits 0** — the graph is a diagnostic and must never abort
a sample (CONSTITUTION principle 7 / rule 8, and consistent with C2 and
C8). On any failure it leaves the un-annotated SVG in place rather than
deleting it or emitting nothing.

Two specific non-fatal conditions, both **warnings**:

- A sentinel that matches no path in the SVG — the node was not drawn.
- A path whose fill matches no sentinel — unexpected, but the picture
  is still valid.

---

## 4. Work items

### 4.1 C12 — `bin/annotate_graph_svg.py` (new)

Container: shares `EXTRACT_BARCODES`'s `neoformit/daff-wf5-scripts`
image. The BandageNG biocontainer has no Python, so C12 cannot run in
it — the same constraint that forced C9 into its own step (rule 14).
Unlike C9 this does **not** need a separate Nextflow process: see §4.2.

Two entry points (subcommands or two scripts, implementer's choice —
one module either way, so the sentinel allocation is defined once and
shared):

- **allocate** — GFA in, `name,color` CSV out.
- **annotate** — SVG + CSV + `assembly_info.txt` in, annotated SVG out.

Parsing notes:

- GFA segment names: `S` records, tab-separated, column 2.
- `graph_path`: `assembly_info.txt` column 8, comma-separated signed
  edge numbers. **Strip the sign** (`-2` and `2` are the same segment,
  traversed in opposite directions) and map number `N` to segment name
  `edge_N`. Treat `*` entries as "no edge" and skip them — they appear
  in real fixture data (`INT-ANIMAL-01` has `*,5,*`).
- SVG editing: prefer a real XML parse over regex where practical, but
  note Qt's output is plain enough that a targeted attribute rewrite is
  defensible if an XML round-trip perturbs the file. Whichever is
  chosen, **the output must still be valid SVG that inlines cleanly
  into the report DOM** — no XML declaration or doctype that would
  break when injected mid-document.

### 4.2 `modules/local/bandage_ng.nf`

Currently emits `${meta.sample_id}.graph.png` from a single
`BandageNG image` call.

Change to the three-step sequence in §3. Because the BandageNG
biocontainer has no Python, the process needs **both** tools available.
Per rule 14, where one process needs multiple biocontainer packages a
**mulled image** should be built (spec 01-pipeline-flow.md §1a). Two
acceptable routes — pick one and record the reasoning:

- **(a)** Build a mulled BandageNG + scripts image, run all three steps
  in one process. Fewer processes, but a new image to pin and maintain.
- **(b)** Split into `BANDAGE_NG` (render, Bandage image) and a small
  `ANNOTATE_GRAPH_SVG` (C12, scripts image), mirroring how
  `MINIPROT_CDS` / `SELECT_GENETIC_CODE` were split for exactly this
  reason (task 38 §2). No new image; one more process and a numbered
  stage-9 sub-step to document.

**(b) is suggested** — it needs no new container, and the precedent is
already established in this pipeline. Note the sentinel CSV must then
be an emitted, staged file crossing the process boundary as a `Path`
(rule 13), not a string.

Update the `stub:` block to `touch` the `.svg` (and the CSV if it is an
emitted output), so `-stub-run` stays green.

`conf/containers.config` needs its `withName: 'BANDAGE_NG'` entry
reviewed against whichever route is taken, and a new entry if (b).

### 4.3 Channel plumbing — `main.nf`

`ch_graph_png` (currently built from `BANDAGE_NG.out.assembly`'s fifth
element) becomes `ch_graph_svg`. `BIN_TARGET` consumes
`BANDAGE_NG.out.assembly` and must keep receiving the same tuple shape
it does now — check the destructuring in the `BIN_TARGET` call and the
`ch_graph_png` map, both of which name the graph element positionally.

### 4.4 `bin/collate.py` (C6)

`--graph-png` → `--graph-svg`; the `diagnostics/` copy target
`graph.png` → `graph.svg`; the `render()` call's `graph_png=` kwarg
renamed. Spec 00-overview.md's bundle layout has already been updated
to say "assembly graph SVG".

### 4.5 Report renderer — `bin/report/report.py`

`graph_png_src` (a base64 `data:` URI via `_file_src`) is replaced by
an inlined SVG string via **`_read_svg`, which already exists** and is
already used for the organelle map — reuse it rather than writing a
second reader (rule 19). Its existing zero-byte/missing handling is
exactly the behaviour wanted here.

Name the context key `graph_svg` for symmetry with
`organelle_map_svg`.

### 4.6 Report template —
`scripts/report/templates/components/assembly.html`

Task 46 placed the graph and the organelle map in a two-column "Genome
annotation" row, each opening a fullscreen modal on click. The graph is
currently an `<img>` with a base64 `src`, duplicated into its modal.

- Swap the `<img>` for the inlined SVG, so tooltips bind (spec
  06a-reports.md's renderer-inputs bullet now groups it with the
  organelle map for this reason).
- **The map's modal already handles the inline-SVG duplication problem
  by *moving* the node into the modal and back rather than copying it**
  — two live copies of one inline SVG collide on internal `id`s and
  `url(#...)` references. The graph SVG now needs the same treatment;
  factor the move/restore behaviour into something both use rather than
  writing it twice.
- Extend the existing info-badge text: say that hovering a node names
  the segment, and state plainly that **node colour carries no binning
  meaning** (§7).

### 4.7 `bin/README.md`

Its component table lists one row per custom-logic script. Add C12
alongside the existing entries, and correct the `report/` and
`render_report.py` rows if the graph rename touches their described
inputs.

### 4.8 Fall-back when the graph is absent

A sample that soft-fails the coverage gate never reaches stage 9. The
Assembly tab's existing "not available for this sample" path must keep
working — confirm against a `fail` fixture, not just an `ok` one.

---

## 5. Tests

### 5.1 Unit — `scripts/tests/test_annotate_graph_svg.py` (new)

C12 is a custom-logic component and targets **100% branch coverage**
(rule 14). Run via `scripts/pytest.sh`. Cover at minimum:

- Sentinel allocation is deterministic — same GFA in twice, identical
  CSV out.
- Allocated colours are unique, and collide with neither `#ffffff` nor
  `#000000`.
- The emitted CSV carries the `name,color` header (§2 item 5 — without
  it BandageNG silently refuses the file).
- GFA with 1 segment, 3 segments, and 0 segments.
- `graph_path` parsing: signed entries (`-2,1,2`) fold to unsigned;
  `*` entries skipped; a segment traversed by two contigs yields both,
  comma-separated; a segment traversed by none yields empty.
- Annotation attaches `data-node`, `data-contigs` and `<title>`, and
  the sentinel fill is **gone** from the output.
- Warning paths: sentinel with no matching path; path with no matching
  sentinel. Both warn, neither raises, exit stays 0.
- Malformed/truncated `assembly_info.txt` and a zero-byte SVG both
  exit 0 leaving the input SVG untouched.

### 5.2 Unit — `scripts/tests/test_report.py` (extend)

The existing `graph_png_src` assertions become `graph_svg`. Add a case
for a graph SVG that is missing and one that is zero-byte, mirroring
the organelle-map cases already there.

### 5.3 Integration — `tests/integration/assertions.sh`

Add a structural assertion that `diagnostics/graph.svg` exists,
is non-zero, and contains at least one `data-node=` attribute for the
assembling samples. This is the only surface that can catch BandageNG
changing its SVG output shape under a future container bump — the unit
tests necessarily work from a captured sample, not a live render.

### 5.4 Fixtures

No new fixture data is needed — `assembly_graph.gfa` and
`assembly_info.txt` already flow through stage 9, and the unit tests
can use small hand-written GFA/`assembly_info.txt` snippets plus a
captured Bandage SVG committed as a test asset.

**A fresh `-profile integration` run is required** to regenerate
`tests/integration/output/`, because the published bundle's filename
changes (`graph.png` → `graph.svg`).

That run also closes a **carry-forward item folded in from
tasks/todo.md** (recorded there by task 46, 2026-09-04):
`tests/integration/output/` is stale (pre-43a) and `report.html` is
zero bytes for all three samples, which is what produced task 46's
"blank sections" review observations. Task 46's own acceptance
criterion 12 deferred the re-run as disproportionate for a
template-only change; this task changes a pipeline stage and needs the
run regardless, so it inherits the obligation. While running it:
confirm `read_qc` and `coverage.recruitment` populate, re-run
`tests/integration/assertions.sh` (expected to go fully green once
`report.html` is regenerated), and re-review the rendered reports
against task 46's presentation changes.

---

## 6. Spec updates — already applied

No further spec work is required; this task implements what the spec
now says. For reference, the design was written into:

- spec/02-stages.md stage 9 — SVG output, the Qt-dump limitation, the
  sentinel round-trip, the header-row requirement.
- spec/02-stages.md §2.2 — new component C12.
- spec/01-pipeline-flow.md — stage 9 node.
- spec/06a-reports.md — Assembly tab (tab 3) and the renderer-inputs
  bullet.
- spec/00-overview.md — bundle layout.

If the implementation diverges from any of the above (in particular if
route (a) is chosen in §4.2, which changes C12's container line),
reconcile the spec in the same change — the constitution requires the
documents to agree.

---

## 7. Explicit non-goal — do not colour nodes by bin

It is tempting to colour graph nodes green/orange/grey to match the
per-contig coverage chart's target/secondary/off-target buckets. **Do
not**, for two independent reasons:

1. **It is not well-defined.** Per §2 item 6, a repeat edge can be
   traversed by both a target and an off-target contig. A plastid
   inverted repeat is exactly this case, and it is not an edge case —
   it is the normal structure of every plastid genome the pipeline
   assembles. Picking one bucket would silently discard the other.
2. **The data does not exist yet at this stage.** Binning is
   `BIN_TARGET`, stage 10. Stage 9 runs before it.

C12 instead exposes every traversing contig in `data-contigs` and
leaves the colouring neutral, so the ambiguity stays visible rather
than being resolved silently (CONSTITUTION principle 7). The binning
classification is read from the per-contig table, which is
unambiguous.

If per-node binning is genuinely wanted later, it is a separate task
with a real design question to answer — probably "colour by the set of
buckets traversed, with a distinct colour for mixed" — and it would
have to move the annotation step after stage 10.

---

## 8. Acceptance criteria

1. `BANDAGE_NG` emits `<sample_id>.graph.svg`; no stage emits
   `graph.png`.
2. Every drawn node in the published SVG carries `data-node`,
   `data-contigs` and a `<title>` child.
3. Hovering a node in the rendered `report.html` shows its segment
   name with **no JavaScript involved** in producing the label.
4. No sentinel colour survives into the published SVG.
5. Re-running the pipeline on unchanged input produces a byte-identical
   SVG.
6. C12 reaches 100% branch coverage under `scripts/pytest.sh`.
7. `flake8` clean on `bin/*.py` and `scripts/tests/*.py`.
8. `nextflow run . -profile stub -stub-run` green.
9. A fresh `-profile integration` run completes;
   `tests/integration/assertions.sh` is fully green, including the new
   `graph.svg` assertion, and `report.html` is non-zero for all three
   samples.
10. A coverage-gate soft-failed sample still renders its Assembly tab
    with the "not available" fallback, no traceback.
11. The graph and the organelle map both open fullscreen without
    breaking each other's inline SVG (shared move/restore, §4.6).
12. Spec and implementation agree — §6 reconciled if §4.2 route (a) is
    taken.

---

## 9. Out of scope

- Per-node bin colouring (§7).
- Click-to-select, pan/zoom, or any graph interaction beyond hover
  labels. The fullscreen modal from task 46 is the zoom affordance.
- Replacing BandageNG, or rendering the graph client-side from the GFA.
  This was considered and rejected in favour of keeping Bandage's
  assembly-aware layout, which is not trivial to reproduce.
- Any change to `BIN_TARGET`'s classification or to the per-contig
  coverage chart.
- task 44_organelle_map.md's real SVG renderer — stage 14 is still a
  stub and stays one here.
