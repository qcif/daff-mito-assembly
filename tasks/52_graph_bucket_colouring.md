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
hover tooltips for the same graph). Independent — either can land
first. They don't share Python: task 53 is report-template/JS only and
does not touch `bin/annotate_graph_svg.py`. They do both edit
`scripts/report/templates/components/assembly.html` (this task: the
graph info badge + colour key; task 53: a tooltip-conversion script
call), so whichever lands second should re-read the other's diff to
that file before starting. This task owns the `<title>` tooltip *text*
(§3.6) — task 53 just re-presents whatever text C12 emits, so the
bucket label added here shows up in task 53's styled tooltip with no
extra work on either side.

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
repeat edge can be shared between the target contig and an off-target
one.

**How this squares with principle 7.** Task 47's objection was to a
*silently resolved* colour — picking one bucket for a segment that
really belongs to several. This task doesn't do that: a segment
traversed by contigs from more than one bucket gets its own explicit
**mixed** colour, and its tooltip names every bucket involved (§3.6).
The ambiguity stays visible. It's just labelled now instead of being
hidden behind uniform grey.

### Reordering is the enabling change — and it also shortens the critical path

`BANDAGE_NG`/C12 only ever needed `assembly_graph.gfa` and
`assembly_info.txt` (from `METAFLYE`, stage 7, or `MEDAKA` when
polishing is on) to do its job. Its position before `BIN_TARGET` is
purely positional, not a real data dependency: `bin_target.py` never
reads the graph SVG. It is handed it only because `main.nf` feeds
`BIN_TARGET` from `ANNOTATE_GRAPH_SVG.out.assembly` (a 5-tuple that
re-emits the assembly files alongside the SVG), and
`modules/local/bin_target.nf` declares a `path(graph_svg)` input
purely to absorb that fifth element.

A side effect of that plumbing today: **every downstream stage
(`BIN_TARGET` → `BLAST_VALIDATE` → `MINIPROT_CDS` → barcodes,
annotation, …) waits on the three-process graph chain**, a diagnostic
that nothing downstream consumes. After this task, `BIN_TARGET` reads
the assembly directly, and the graph chain hangs off `BIN_TARGET` as a
side branch that runs in parallel with stages 11–14. The graph's
position in the DAG moves, but it stops being on the critical path.

### Why "after `BIN_TARGET`", not "at the end of the run"

The idea was first raised as moving the graph stage "near the end of
the run". The binning classification is finalised at stage 10, and
nothing downstream (BLAST validation, CDS annotation, barcode
extraction) changes *which bucket a contig is in* — those stages use
the binning decision, they don't revise it. Running the graph any
later than immediately after `BIN_TARGET` would add nothing to the
colouring and would chain a diagnostic behind stages it doesn't depend
on. Recommendation: swap stages 9 and 10 only.

---

## 2. What the current code already gives this task

Verified against the current tree:

- `bin/report/report.py`'s `_contig_bucket(classification)` maps
  `bin_metadata.json`'s per-contig `classification` field to one of
  three buckets: `target_candidate` → `target`, `secondary_target` →
  `secondary`, anything else (`off_target`, `sibling_organelle`,
  `None`) → `off-target`. `_BUCKET_COLOURS` pins the palette:
  `target` `#2ca02c` (green), `secondary` `#ff7f0e` (orange),
  `off-target` `#7f7f7f` (grey) — matplotlib/Plotly's tab10 values.
- **The palette is already duplicated.** The coverage chart in
  `scripts/report/templates/components/assembly.html` does not read
  `coverage_chart.colors` (which `_coverage_chart_data()` computes and
  ships); it hardcodes the same three hexes again in a JS
  `bucketOrder` array. Adding the graph as a third consumer without
  fixing this would mean three copies (rule 19). §3.5 fixes it.
- `bin_target.py` writes `bin_metadata.json` with a `contigs` list
  holding **every** contig — `primaries + secondaries`, where
  `secondaries` is every row not classified `target_candidate`
  (`select_primary()`), with `contig_id` and `classification` on each.
  `contig_id` is the Flye contig name (`contig_N`), i.e. the same
  namespace as `assembly_info.txt` column 1 and therefore as
  `data-contigs`. It also writes `contigs_selected` (the target
  contig ids).
- `bin/annotate_graph_svg.py`'s `parse_graph_path_column()` already
  returns `{segment_name: {contig_id, ...}}` — the join key against
  those classifications.
- `bin/organelle_map.py`'s `load_bin_metadata()` (task 44) is the
  existing "absent/empty/malformed → `{}`" loader for the same file —
  the degradation shape to copy.
- `modules/local/bin_target.nf` already emits `tuple val(meta),
  path("bin_metadata.json"), emit: metadata`, joinable by `meta` — no
  new Nextflow output needed. Its stub `touch`es an empty
  `bin_metadata.json`, which exercises the degrade path under
  `-stub-run` for free.
- `bin/annotation_gff.py` (shared between C10/C11/`organelle_map.py`,
  imported as a plain top-level module from `bin/`) is the precedent
  for a small shared module in `bin/` with no Nextflow process of its
  own. Scripts in `bin/` get `bin/` as `sys.path[0]`, and every
  `scripts/tests/test_*.py` already inserts `BIN_DIR` on `sys.path`,
  so a new top-level `bin/contig_bucket.py` is importable from
  `annotate_graph_svg.py`, `bin/report/report.py` (run via
  `bin/render_report.py`) and the tests without any packaging change
  or container rebuild.

---

## 3. Design

### 3.1 Stage renumbering

Swap stage 9 and stage 10 throughout the spec and code comments:
`BIN_TARGET` becomes stage 9, `BANDAGE_NG` (+ C12) becomes stage 10.
No other stage number changes. Update:

- spec/01-pipeline-flow.md's flow diagram (the `[9] BANDAGE_NG` /
  `[10] BIN_TARGET` lines and their arrows — the graph now branches
  off after binning).
- spec/02-stages.md:
  - the stage table rows for 9/10 (swap order and numbers);
  - the `BANDAGE_NG` row's closing sentences, which currently say
    **"Bucket colouring of graph nodes is deliberately not
    specified"** and "Binning is not known at this stage anyway
    (`BIN_TARGET` is stage 10)" — rewrite to describe §3.4's
    colouring rule, including the mixed colour and why it satisfies
    principle 7 (§1);
  - the C12 row in §2.2: its "Stage" link, and its **"Never invents a
    per-node bucket … the colouring stays neutral"** sentence —
    rewrite the same way. Its inputs gain `bin_metadata.json`. The
    three-process split is unchanged; only its position in the DAG
    moves;
  - the C3 row's stage reference.
- spec/06a-reports.md's Assembly tab row (§3): the "node colour
  carries no binning meaning … the binning classification is read
  from the per-contig table, not from the picture" sentence is no
  longer true. Rewrite it to say what colour means now, and that mixed
  is shown explicitly rather than resolved. The "[§2 stage 9]" link
  text becomes stage 10.
- `bin/annotate_graph_svg.py`'s module docstring: "Stage 9" → "Stage
  10", and the **"Never invents a per-node bin colour"** paragraph →
  describe the bucket/mixed/unknown rule and the `--bin-metadata`
  input.
- Header comments in `modules/local/bandage_ng.nf`,
  `modules/local/allocate_graph_sentinels.nf`,
  `modules/local/annotate_graph_svg.nf` and
  `modules/local/bin_target.nf`, plus the `// Stage 9` / `// Stage
  10` / `// Stages 7–9` comments in `main.nf`.
- `bin/README.md`'s component table (C3: `10 BIN_TARGET` → `9
  BIN_TARGET`; C12: `9 ALLOCATE_GRAPH_SENTINELS, ANNOTATE_GRAPH_SVG` →
  `10 …`).
- `tests/integration/assertions.sh`'s progressive-uncomment-plan
  comment block at the top (the `BANDAGE_NG (real renderer)` /
  `BIN_TARGET (real C3)` rows) — and move the `graph.svg` assertion
  block below the `BIN_TARGET` block so the file reads in stage order.

Check for any other "stage 9"/"stage 10" references before finishing
(`grep -rn -i "stage 9\|stage 10" spec bin modules main.nf tests`).
Leave completed task files in `tasks/completed/` alone — they describe
the pipeline as it was.

### 3.2 `main.nf` reordering

Sketch (not literal code):

```
ch_assembly = params.polish ? MEDAKA(...).assembly : METAFLYE.out.assembly

// Stage 9: bin contigs — reads the assembly directly again
BIN_TARGET(ch_assembly, ch_organelle_refs)

// Stage 10: graph diagnostic — off the critical path, joined to
// BIN_TARGET's classification only at the final annotate step
ALLOCATE_GRAPH_SENTINELS(ch_assembly)           // unchanged
BANDAGE_NG(ALLOCATE_GRAPH_SENTINELS.out.sentinels)  // unchanged
ANNOTATE_GRAPH_SVG(
    BANDAGE_NG.out.rendered.join(BIN_TARGET.out.metadata, by: 0))

// Stage 11 onwards: unchanged, no longer waits on the graph chain
BLAST_VALIDATE(BIN_TARGET.out.binned)
```

Note `ALLOCATE_GRAPH_SENTINELS` and `BANDAGE_NG` don't need binning,
so they still start as soon as the assembly exists and run
concurrently with `BIN_TARGET`. Only `ANNOTATE_GRAPH_SVG` waits on the
join.

**Narrow `ANNOTATE_GRAPH_SVG`'s output.** It currently re-emits
`(meta, assembly, gfa, info, graph.svg)` only so `BIN_TARGET` could
consume it. Once `BIN_TARGET` no longer does, the only consumer is
`ch_graph_svg`, which immediately `.map`s it down to `(meta,
graph_svg)`. Change the output to `tuple val(meta),
path("${meta.sample_id}.graph.svg"), emit: graph` and drop the `.map`
in `main.nf`. `ch_graph_svg`'s shape for `COLLATE` stays the same.

**Join cardinality.** The join is by `meta`, and both sides emit
exactly once per assembling sample, so no `remainder:` is needed. A
sample that fails `BIN_TARGET` outright (process error, not an empty
bin) will now also lose its graph. Today it would have lost everything
downstream anyway, including `COLLATE`'s `ch_ok_inputs` join, so
nothing new is lost. Note this in the module comment rather than
working around it.

### 3.3 Module changes

- `modules/local/annotate_graph_svg.nf`: input tuple becomes
  `(meta, assembly, gfa, info, sentinels_csv, raw_svg,
  bin_metadata_json)` — or narrow the pass-through from `BANDAGE_NG`
  too, if that reads cleaner; implementer's call. Pass
  `--bin-metadata ${bin_metadata_json}` to `annotate_graph_svg.py
  annotate`. Stub unchanged.
- `modules/local/bin_target.nf`: drop the unused `path(graph_svg)`
  from the input tuple.
- `modules/local/bandage_ng.nf`, `allocate_graph_sentinels.nf`: comment
  changes only (§3.1).

### 3.4 `bin/annotate_graph_svg.py` — bucket colouring

New `annotate` argument: `--bin-metadata PATH` (optional, default
`None`). Loaded with the same "missing / empty / malformed JSON →
`{}`" shape as `organelle_map.load_bin_metadata()`. Don't import
`organelle_map` for it — one small loader in `contig_bucket.py`
(§3.5), used by both, is the rule-19 answer if the implementer wants
to deduplicate; a local copy is acceptable if it isn't shared.
Build `{contig_id: classification}` from its `contigs` list (skipping
entries with no `contig_id`).

New pure function, used by `annotate_svg()` in place of the current
always-`DISPLAY_COLOUR` fill:

```
segment_bucket(contig_ids, contig_classifications) -> str | None
    buckets = { contig_bucket(contig_classifications[c])
                for c in contig_ids if c in contig_classifications }
    none of the traversing contigs is classified   -> None
    more than one distinct bucket                  -> 'mixed'
    exactly one                                    -> that bucket
```

Colour choice per node:

| Segment state | Fill | Tooltip suffix (§3.6) |
|---|---|---|
| only `target` contigs | `BUCKET_COLOURS['target']` `#2ca02c` | `(target)` |
| only `secondary` | `BUCKET_COLOURS['secondary']` `#ff7f0e` | `(secondary)` |
| only `off-target` | `BUCKET_COLOURS['off-target']` `#7f7f7f` | `(off-target)` |
| ≥2 distinct buckets | `MIXED_COLOUR` — `#9467bd` (tab10 purple) | `(mixed: target, off-target)` — buckets in fixed target/secondary/off-target order |
| no classified traversing contig, or no `bin_metadata.json` | existing `DISPLAY_COLOUR` `#c8c8c8` | none — tooltip unchanged from task 47 |

Notes:

- **Why `#9467bd` for mixed:** it continues the tab10 palette the other
  three already come from, it isn't used anywhere else in the report,
  and it reads clearly against white and against both greys. Declare
  it as a constant next to `BUCKET_COLOURS` in `contig_bucket.py`, not
  inside `BUCKET_COLOURS` (so `report.py`'s coverage chart, which only
  ever has single-bucket contigs, can't pick it up by accident).
- **Two greys.** Off-target (`#7f7f7f`, mid grey) and unknown
  (`#c8c8c8`, light grey — task 47's neutral) are deliberately
  different: "the pipeline looked and this is not the target" is not
  the same statement as "the pipeline has no classification for this
  segment". The colour key (§3.7) must list both. In a normal run,
  unknown only appears for a segment no contig's `graph_path`
  traverses, or when `bin_metadata.json` is missing — the report
  should make neither look like an off-target call.
- A partially classified segment (some traversing contigs classified,
  some not) is coloured from the classified ones only. It's rare
  enough (would need `bin_metadata.json` to omit a contig Flye
  reported) that a separate state isn't worth it — but the
  unclassified contig ids are still listed in `data-contigs` and the
  tooltip, so nothing is hidden.
- Add a `data-bucket` attribute (`target` / `secondary` /
  `off-target` / `mixed`, omitted when unknown) alongside
  `data-node`/`data-contigs`. It gives the integration assertion
  (§5.4) and any later report styling a stable hook that doesn't
  depend on matching hex strings.
- **Never fails the sample.** `annotate` already always exits 0 and
  falls back to the unannotated SVG on any exception (CONSTITUTION
  principle 8, rule 8); the new input must stay inside that envelope.
  A missing/malformed `bin_metadata.json` is not an exception at all —
  it is the "unknown" row above, i.e. today's output exactly.

`_recoloured_open_tag()` takes the chosen fill colour and bucket label
as parameters instead of hardcoding `DISPLAY_COLOUR`. The sentinel
matching/round-trip logic is untouched.

### 3.5 `bin/contig_bucket.py` (new, shared)

Single definition (rule 19) of:

- `contig_bucket(classification) -> str` — moved verbatim from
  `report.py`'s `_contig_bucket`.
- `BUCKET_COLOURS` — moved from `report.py`'s `_BUCKET_COLOURS`.
- `BUCKET_ORDER` (`target`, `secondary`, `off-target`) and display
  labels (`Target` / `Secondary` / `Off-target`) — currently only in
  the template's JS.
- `MIXED_COLOUR` (§3.4). The unknown colour stays
  `annotate_graph_svg.DISPLAY_COLOUR` — it's C12's own neutral, not a
  bucket.
- Optionally `load_bin_metadata()` (§3.4).

It has to be a top-level `bin/` module, not something under
`bin/report/`: importing `report.*` from `annotate_graph_svg.py` would
execute the report package's imports (`jinja2` etc.), which the
`wf5-scripts` image running `ANNOTATE_GRAPH_SVG` has no reason to
carry.

Consumers:

- `bin/report/report.py`: import `contig_bucket`/`BUCKET_COLOURS`
  instead of defining them; delete `_contig_bucket`/`_BUCKET_COLOURS`.
- `scripts/report/templates/components/assembly.html` coverage chart:
  replace the hardcoded `bucketOrder` JS array with data from
  `report.py` (e.g. `assembly_view.bucket_legend` — a list of `{key,
  name, color}` built from `BUCKET_ORDER`/labels/`BUCKET_COLOURS`,
  emitted with `| tojson`). This removes the third copy of the
  palette. It must not change what the chart renders.

### 3.6 Tooltip text

C12 owns the `<title>` text (task 47 §4.4). Extend it with the bucket
suffix from the §3.4 table so the colour explains itself on hover —
most importantly in the **standalone** `diagnostics/graph.svg`, which
has no legend or info badge around it:

```
edge_3 — contig_1 (target)
edge_7 — contig_1,contig_4 (mixed: target, off-target)
edge_9                       # no traversing contig: unchanged
```

Task 53 re-presents this text as a Bootstrap tooltip in the report,
so no coordination is needed beyond task 53 not re-deriving the text
(its §3.3 already says so).

### 3.7 Report template: info badge + colour key

`scripts/report/templates/components/assembly.html`:

- Rewrite the "Assembly graph" `info_badge` text. It currently ends
  "node colour is not meaningful and does not indicate
  target/secondary/off-target classification — see the contig table
  above for that". Replace that with a plain-language description of
  the five states. Use the same bucket words as the coverage chart
  legend and the contig table's Classification column, so a reader
  doesn't have to reconcile two descriptions of one palette.
- Add a small static colour key directly under the graph thumbnail:
  a plain HTML swatch list (five inline swatches + labels, rendered
  from the same `bucket_legend` data as §3.5 plus the mixed/unknown
  entries), not a Plotly legend — the graph is an inline SVG, not a
  Plotly figure, and a swatch list works identically in the thumbnail
  and needs no JS. Only list states that actually occur in this
  sample's graph if that's cheap to compute (`data-bucket` values,
  §3.4), otherwise list all five.
- The key must sit **outside** the `#assembly-graph` div: the
  modal move/restore helper (task 47 §4.6) moves that div's children
  into the fullscreen modal, and the key shouldn't disappear from the
  page when the modal opens. Also render it inside the graph modal's
  body, below the moved SVG, so the fullscreen view has it too.

---

## 4. Work items

1. `bin/contig_bucket.py` (new) per §3.5; switch `bin/report/report.py`
   to it; replace the coverage chart's hardcoded `bucketOrder` with
   `report.py`-supplied legend data.
2. `bin/annotate_graph_svg.py`: `--bin-metadata`, `segment_bucket()`,
   per-node fill + `data-bucket` + tooltip suffix per §3.4/§3.6;
   docstring rewrite per §3.1.
3. `main.nf` reorder + `ANNOTATE_GRAPH_SVG` output narrowing (§3.2);
   module input changes (§3.3).
4. `assembly.html` info badge + colour key (§3.7).
5. Spec, `bin/README.md`, module/`main.nf` comment and
   `assertions.sh` comment reconciliation (§3.1).
6. Stub run, unit tests, fresh integration run, rendered-report
   hand-check (§5).

---

## 5. Tests

### 5.1 Unit — `scripts/tests/test_contig_bucket.py` (new)

100% branch coverage (rule 14):

- `contig_bucket()` for every `classification` value `bin_target.py`
  actually emits (`target_candidate`, `secondary_target`,
  `off_target`, `sibling_organelle`), plus `None` and an unrecognised
  string (both → `off-target`).
- `BUCKET_COLOURS` keys == `BUCKET_ORDER`; `MIXED_COLOUR` not among
  `BUCKET_COLOURS` values and ≠ `annotate_graph_svg.DISPLAY_COLOUR`.
- `load_bin_metadata()` (if it lives here): valid file, missing path,
  `None`, empty file, malformed JSON.

### 5.2 Unit — `scripts/tests/test_annotate_graph_svg.py` (extend)

- `segment_bucket()`: each single-bucket case; two-way and three-way
  mixed; a contig in `contig_ids` absent from the classifications
  (ignored, not an error); all contigs unclassified → `None`; empty
  `contig_ids` → `None`.
- `annotate_svg()`: node fill is the bucket colour; a mixed node gets
  `MIXED_COLOUR` (not a bucket colour, not the neutral one);
  `data-bucket` present with the right value, absent for unknown;
  `<title>` suffix format per §3.6, including the fixed bucket order in
  the mixed case; unknown node's `<title>` is byte-identical to task
  47's.
- `run_annotate()` / `main()`: `--bin-metadata` omitted, pointing at a
  missing file, empty file and malformed JSON — all produce today's
  all-neutral output and exit 0. Existing task 47 tests must keep
  passing, with only the changes the new signature forces.

### 5.3 Unit — `scripts/tests/test_report.py` (update)

- Point the `test_contig_bucket_*` cases at `bin/contig_bucket.py`, or
  move them into §5.1 and delete them here — don't keep two copies.
- Coverage chart assertions keep passing unchanged: this task must not
  change what the chart shows.
- New: the rendered Assembly tab contains the colour key, and the
  graph info badge no longer contains "not meaningful".

### 5.4 Integration — `tests/integration/assertions.sh`

Extend the existing `graph.svg` block (task 47 §5.3) with an invariant
that needs no hardcoded edge numbers and no `graph_path` parsing in
bash (`data-contigs` already carries the traversal):

- For each id in `bin_target/bin_metadata.json`'s `contigs_selected`,
  every `data-node` group whose `data-contigs` list contains it has
  `data-bucket` of `target` or `mixed`.
- At least one node in the SVG has `data-bucket="target"` or
  `"mixed"` — i.e. the target is visible in the picture at all.
- Every `data-node` group's `fill` is one of the five §3.4 colours —
  no sentinel fill leaks through now that the recolour path branches.

Keep the structural checks from task 47 as they are.

### 5.5 Stub and fixtures

- `nextflow run . -profile stub -stub-run` must pass with the reordered
  DAG (empty stub `bin_metadata.json` → degrade path).
- No new fixture data: `INT-ANIMAL-01`, `INT-PLANT-01-pt` and
  `INT-PLANT-01-mt` already produce real `bin_metadata.json`.
- A fresh `-profile integration` run (clean `./work/` — the DAG
  changed) must leave `tests/integration/assertions.sh` fully green.
  Render a report and hand-check that the Assembly tab graph shows
  coloured nodes that agree with the contig table and coverage chart,
  and that the key reads correctly. Give the report path to the user
  (CLAUDE.md).

---

## 6. Explicit non-goals

- Changing what `BIN_TARGET`'s classification means or its thresholds
  (`params.bin_target_thresholds`) — this task only *displays* the
  existing decision on a second picture.
- Resolving `mixed` segments to a single bucket by any heuristic
  (length share, coverage, …) — that is exactly the silent resolution
  principle 7 forbids.
- Distinguishing `sibling_organelle` from `off_target` in the graph.
  Both map to `off-target` today in the chart and table, and the graph
  must agree with them. If that distinction is wanted, it's a
  `contig_bucket.py` change that would apply to all three views at
  once — a separate task.
- A legend baked into the standalone `graph.svg` — the tooltip suffix
  (§3.6) is the standalone file's self-description.
- Interaction beyond the static key (click-to-filter, bucket toggles)
  — same boundary task 47 §9 drew.
- Colouring the organelle map (task 44_organelle_map.md) by bucket —
  it colours by feature kind (CDS/tRNA/rRNA), an unrelated axis.

---

## 7. Acceptance criteria

1. `BIN_TARGET` consumes the assembly directly; the graph chain
   branches off after it and no longer blocks `BLAST_VALIDATE` or
   anything downstream. Stage numbers (9 = `BIN_TARGET`, 10 =
   `BANDAGE_NG` + C12) are consistent across spec, `main.nf`, module
   comments, `bin/README.md` and `assertions.sh` comments.
2. Graph node fills follow the §3.4 table exactly: target `#2ca02c`,
   secondary `#ff7f0e`, off-target `#7f7f7f`, mixed `#9467bd`,
   unknown `#c8c8c8`. Each classified node carries `data-bucket`, and
   its `<title>` carries the matching suffix.
3. `bin_metadata.json` missing, empty or malformed does not fail the
   sample (exit 0); every node falls back to task 47's neutral output.
4. `bin/contig_bucket.py` is the only definition of the
   classification → bucket → colour mapping. `report.py`, the coverage
   chart template and `annotate_graph_svg.py` all read from it; no
   bucket hex remains hardcoded in `assembly.html`.
5. New/changed `bin/*.py` at 100% branch coverage under
   `scripts/pytest.sh`; `flake8` clean.
6. `nextflow run . -profile stub -stub-run` green with the reordered
   DAG.
7. A fresh `-profile integration` run completes;
   `tests/integration/assertions.sh` is fully green, including the
   §5.4 checks.
8. The Assembly tab's info badge and colour key describe what node
   colour now means. Nothing in spec/02-stages.md, spec/06a-reports.md,
   the C12 docstring or the template still claims graph colour is
   meaningless or deliberately neutral.
