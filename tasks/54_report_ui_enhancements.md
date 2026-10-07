# Task 54 — Per-sample report UI enhancements

**Phase:** P4a (spec 06-phases.md).
**Goal:** Act on a design review of the per-sample `report.html`:
tighten the layout, replace number-heavy tables with small charts where
a picture reads faster, fix two rendering defects (oversized assembly
graph, clipped genome-map legend), and add two small features (assembly
FASTA copy/download, barcode-locus labels on the genome map).

**Reviewed report:** `reference-material/task-54-report-int-animal-01.html`
(an `INT-ANIMAL-01` render). Use it as the "before" for every visual check
below.

**Depends on:** task 52_graph_bucket_colouring.md. Its changes to
`scripts/report/templates/components/assembly.html` and
`bin/report/report.py` must be committed before this task starts.

**Related tasks:** task 53_graph_bootstrap_tooltips.md has landed
(commit d1ec7bd). It added `scripts/report/static/js/svg-title-tooltips.js`
and the graph's `promoteSvgTitleTooltips` /
`hideSvgTooltipsAroundModal` calls in `assembly.html`. Those calls move
together with the graph div and its modal (§5.2). §5.5 applies the same
helpers to the genome map.

---

## 1. Overview

The per-sample report (spec §6a) is the main thing a biosecurity
officer reads. A review of a real `animal_mt` render found these
problems:

- **Overview** is laid out as one long column. *Key findings* and
  *Warnings* belong side by side. The *Key findings* table shows
  severity by tinting whole rows, which is hard to scan and drowns out
  the values. Read QC stacks a table above a wide horizontal chart.
- **Validation** shows recruitment and coverage as a long table of
  large numbers. The reader really asks three questions: what fraction
  of reads did we keep, what fraction of bases, and where does the
  coverage estimate sit relative to the gate floors? Each question is
  easier to answer as a small bar chart. Subsampling detail is rarely
  relevant, but it takes up four rows every time. Genetic codes show as
  bare NCBI table numbers ("5"), which mean nothing to most readers.
- **Assembly**: the Bandage graph sits under *Genome annotation*, but
  it shows how the assembly is put together and has nothing to do with
  annotation. It also renders very large and overflows its column. The
  genome map's legend footnote is clipped at the bottom of the SVG. The
  map does not show which of its many genes are the barcode loci, even
  though those are the reason the pipeline exists. The assembled
  organelle FASTA cannot be reached from the report at all.
- **Barcodes** always shows a "Panel barcodes not recovered" section,
  even when it is empty.

None of this changes a result, a status or a threshold. It is
presentation only. §5.4 gives `ORGANELLE_MAP` the locus-panel config
file so it can label barcode loci, but the stage's channel inputs do
not change.

**Constraints that still apply to every change below:**

- The report stays a self-contained single file (spec §6a.5). No new
  external assets. Icons are inline SVG or Unicode, following the
  inline Font Awesome SVG already used in
  `scripts/report/templates/components/sequence-display-modal.html`.
- Negative clarity (CONSTITUTION principle 7). `low_coverage` must
  never look like a failure, and `no_assembly` and `no_barcode` must
  never look alike. This applies to every new badge and colour.
- Auditability (CONSTITUTION rule 18). When a chart replaces a table
  value, the exact figure must stay readable **without hovering**. Show
  it as a bar text label or a caption, because a printed or saved
  report has no hover. The same rule covers genetic codes: the
  human-readable name replaces the raw ID on screen, and the NCBI table
  number stays visible in a tooltip.
- Severity is still computed in Python (`bin/report/report.py`, which
  sets `class` on each finding), never in a template (task 46 §7).

---

## 2. What the current code gives this task

These points were checked against the current tree:

- **Key findings**: `scripts/report/templates/components/key-findings.html`
  renders `key_findings` as `<tr class="table-{{ finding.class }}">`.
  Severities in use are `success`, `warning`, `danger`, `info` (the
  `low_coverage` outcome) and `secondary` (assembly size, unknowns).
  `scripts/tests/test_report.py` asserts `table-info` on the
  low-coverage row (around the `id="key-findings"` assertions), so that
  test must change with the markup.
- **Overview layout**: `components/overview.html` includes
  key findings, then an inline warnings block, then `read-qc.html`.
- **Read QC**: `components/read-qc.html` contains the table and a
  horizontal Plotly stacked bar (`#filter-yield-chart`, 500×120 px).
- **Validation**: `components/validation.html`. The recruitment/gate
  table sits in a `col-md-6` beside *Annotation cross-checks*. Data
  comes from `validation_view` (`report.py`): `recruitment`,
  `raw_bases`, `gate.*`, `flye_depth.coverages` (one value per target
  contig) and `estimate.*` (subsampling).
- **Genetic codes** appear in three places: Validation
  (`genetic_code_annotate`, `genetic_code_cds`) and the Barcodes table
  (`locus.genetic_code`). The NCBI tables configured per target are
  listed in spec/03-organelles.md (§3.1 table: 11, 1, 2/5) and in
  `params.genetic_code_tables`.
- **Assembly**: `components/assembly.html`. *Assembly statistics* is a
  `col-sm-4` table plus a `col-sm-8` per-contig coverage chart. The
  graph (`#assembly-graph`, its colour key, `#graphModal`) is in the
  `col-md-4` of *Genome annotation*, beside the `col-md-8` *Genome map*
  `<h4>`.
- **Assembly FASTA is not in the report.** `COLLATE` already copies
  `target_fasta` to `organelle_assembly.fasta` (`bin/collate.py`), but
  `render()` in `bin/report/report.py` takes no FASTA argument except
  `barcodes_fasta`. `_file_src()` already base64-encodes a file into a
  `data:` URI. Reuse it.
- **Genome map clipping**: in `bin/organelle_map.py`'s `render_svg()`,
  the legend starts at `height - 4*LEGEND_ROW_HEIGHT + 10`. It has three
  kind rows and then the strand footnote at
  `+ 3*LEGEND_ROW_HEIGHT + 14`, which puts the footnote baseline at
  `height + 6`, below the viewBox. That is the cut-off text.
- **Genome map gene names**: `ORGANELLE_MAP` already draws every CDS
  with its gene symbol (`locus['gene']`, shown in the tooltip). The
  barcode panel is `assets/loci.json` (`params.locus_panel`, keyed by
  `assembly_target`). `bin/validate_barcodes.py` selects panel loci
  with exactly this name match: panel symbol `.upper()` against the
  upper-cased CDS gene symbol. It receives the panel inline as
  `${file(params.locus_panel)}` in `modules/local/extract_barcodes.nf`.
- **Barcodes**: `components/barcodes.html` renders the
  `barcodes_view.dropped` section unconditionally, with an "Every panel
  locus was recovered." empty row.

---

## 3. Overview tab

### 3.1 Key findings + warnings: two columns

- Put *Key findings* and *Warnings* in one Bootstrap row, at
  `col-md-8` / `col-md-4`. Below `md` they stack, Key findings first.
- *Run context* (inputs, wall time, provenance) stays above both, as
  spec §6a.2 tab 1(a) requires.

### 3.2 Severity badges instead of row colours

- Remove the `table-{class}` row tint. Each value cell gets a small
  round badge **to the left of the text**, coloured and iconed by
  severity:

  | `class` | Badge | Icon |
  |---|---|---|
  | `success` | green | tick |
  | `warning` | orange | `!` |
  | `danger` | red | `X` |
  | `info` | blue | `i` |
  | `secondary` | grey | none / neutral dash |

  The review asked for only the first three. `info` and `secondary`
  are added because the data already uses them. `info` carries the
  `low_coverage` outcome, which must not get a warn or danger icon
  (principle 7).
- Implement the badge as a new macro (e.g.
  `scripts/report/templates/macros/severity-badge.html`) so the
  run-level report can reuse it later. Give each badge a
  `title`/`aria-label` naming the severity in words, so meaning does
  not depend on colour alone.
- Keep `id="key-findings"` on the table. Tests anchor on it.

### 3.3 Top BLAST hit value in a smaller font

The *Top BLAST hit* value (`"98.28% — <long stitle>"`) wraps over
several lines. Render that row's value in `small` text. To keep
severity logic out of the template, let `_top_blast_hit_finding()` set
a presentation hint (e.g. `'compact': True`) and have the template
apply the class from that hint. Do not test the label string in the
template.

### 3.4 Sequencing data quality

- Table and chart sit side by side, in two equal columns (`col-md-6`
  each), stacking below `md`.
- Rename the table header "Clean reads (post CHOPPER + FILTLONG)" to
  **"Clean reads"**. The filter-yield sentence under the table still
  names CHOPPER + FILTLONG, so no information is lost.
- Change the filter-yield chart to **vertical** stacked bars:
  categories Reads / Bases on x, retained (green) and discarded (grey)
  percentages on y, range 0–100. Show the retained percentage as a bar
  text label (rule 18). Let the chart fill its column instead of using
  the fixed 500 px width.

---

## 4. Validation tab

### 4.1 Gate decision as a headline badge

Move *Gate decision* out of the table into a large badge component at
the top of the *Recruitment & coverage gate* section (e.g. a
`fs-5`-sized badge, or a small bordered card holding the badge and a
one-line meaning). Keep the existing colour mapping: `ok` success,
`low_coverage` warning, `fail` danger. Keep the info badge text. The
`low_coverage` warn-floor alert and the `total_recruited` alert stay
directly under it. `tests/integration/assertions.sh` greps the
warn-floor alert text on `INT-PLANT-01-mt`, so do not reword it.

### 4.2 Three small charts replace six table rows

Add a row of three small vertical Plotly bar charts. Each chart shows a
**subset against its whole**: grey for the whole, green for the subset.

| Chart | Grey (whole) | Green (subset) | Replaces rows |
|---|---|---|---|
| Reads | reads aligned | reads recruited | "Reads aligned / recruited" |
| Bases | total raw bases | total recruited bases | "Total raw bases", "Total recruited bases" |
| Coverage | estimated coverage, recruited reads | estimated coverage, assembled organelle (Flye) | both "Estimated coverage" rows |

The review says "stacked". Do **not** literally stack the two values:
the green figure is a subset of the grey one, so a stack would
double-count. Draw them as nested bars instead: a wide grey bar with a
narrower green bar overlaid on it (Plotly `barmode: 'overlay'` with
different `width`s, or an equivalent). This reads as "green out of
grey", which is what the review means. The coverage pair is not a strict
subset, because assembled depth can exceed the recruited estimate. A
narrower green bar stays visible in either case.

- **Coverage floors**: on the coverage chart, draw the hard-min, warn
  and max floors (`gate.hard_min_required`, `gate.warn_threshold`,
  `gate.max_allowed`) as labelled horizontal lines (Plotly `shapes` +
  `annotations`). Use colours that match the gate badge (danger /
  warning / neutral). The y-axis must include all three floors, so the
  reader can see where the estimate sits against each.
- **Multiple target contigs** (`flye_depth.coverages` has more than one
  value, e.g. multi-contig `plant_mt`): draw one narrow green bar per
  target contig in the coverage chart, labelled by contig. Do not
  aggregate them into a derived number the pipeline never computed
  (rule 18).
- **Missing data**: a chart whose inputs are absent (e.g. no
  `recruitment` block, no `flye_depth` on a soft-failed sample) shows
  only the bars it has. If none of its inputs exist, render a muted
  "not available" line instead of the chart. Never draw a zero-height
  bar for a missing value, because that is the confident-negative
  rendering principle 7 forbids.
- Show the exact values as bar text labels, formatted the way the table
  formats them now (thousands separators, `×` suffix). Keep each removed
  row's info-badge explanation, moved to the chart's title or caption.
- The floors row ("Fail / warn / max coverage floors") and the plant
  sibling-split rows stay in the table.

### 4.3 Subsampling in a modal

- Remove the four subsampling rows (fraction, seed, pre/post coverage)
  from the main table. Put them in a modal (`#subsamplingModal`) with
  the same info badges.
- In the main table, a single row reads **Subsampled: No**, or
  **Subsampled: Yes** followed by a "Subsampling" button that opens
  the modal. The button and the modal are rendered only when
  `estimate.subsampled` is true. Keeping the row means a reader can
  still confirm that subsampling *did not* happen (rule 18).

### 4.4 Human-readable genetic codes

- Add a constant mapping from NCBI translation-table ID to a name in
  `bin/report/report.py` (or `bin/report/config.py` if that is where
  report constants live). Cover at least every table in
  `params.genetic_code_tables` (1, 2, 5, 11) plus 4 and 9. Use the NCBI
  names in short form, e.g. 1 "Standard", 2 "Vertebrate mitochondrial",
  5 "Invertebrate mitochondrial", 11 "Bacterial, archaeal and plant
  plastid".
- Expose it as a Jinja filter (e.g. `genetic_code_name`) and apply it
  in Validation (both rows) and the Barcodes table.
- Render the name, with a tooltip that reads "NCBI translation table
  N". An ID not in the mapping renders as "NCBI table N" and is never
  blank. `None` still renders as `-`.

---

## 5. Assembly tab

### 5.1 Assembly FASTA: copy and download

- Add an `assembly_fasta: Optional[Path]` argument to `render()` /
  `build_context()`. Encode it with `_file_src()` as
  `assembly_fasta_src`, and treat a missing or zero-byte file as `None`,
  the same way `annotation_gff` is handled. Pass `args.target_fasta`
  from `bin/collate.py`, and add the matching flag to
  `bin/render_report.py`.
- In *Assembly statistics*, add an "Assembly FASTA" line with two icon
  buttons: **copy to clipboard** (reuse the copy icon and feedback
  behaviour from `sequence-display-modal.html`) and **download**
  (`download="organelle_assembly.fasta"`). Give both a tooltip and an
  `aria-label`. The tooltip must state that this is the **target**
  contig(s) only, not the whole metaFlye assembly.
- Embed the sequence once. Copy-to-clipboard decodes the existing
  `data:` URI in JS (`atob`) instead of inlining the sequence a second
  time. A plant mitogenome can be several hundred kb, so doubling it
  would be wasteful.
- No `COLLATE` process or `main.nf` change is needed: `target_fasta` is
  already a `COLLATE` input.

### 5.2 Move the assembly graph into *Assembly statistics*

- *Assembly statistics* becomes three columns: stats table
  `col-sm-4`, per-contig coverage chart `col-sm-5`, assembly graph
  `col-sm-3`. Fit the column widths to the reviewed report if 4/5/3
  turns out cramped, but the graph column stays `col-sm-3` as
  requested.
- Move the graph `<h4>`, its info badge, the colour key, the
  click-to-zoom wrapper, `#graphModal` and its
  `bindSvgModalMoveRestore` call together. The colour key must stay
  *outside* `#assembly-graph`; the existing template comment explains
  why.
- **Fix the overflow.** The inlined Bandage SVG renders far larger than
  its column. Find the cause (likely fixed `width`/`height` attributes
  on the root `<svg>` from BandageNG, or no size constraint on the
  container) and fix it in the report: CSS in
  `scripts/report/static/css/main.css` scoped to
  `#assembly-graph svg` (`max-width: 100%; height: auto;`), or by
  removing/overriding the root size attributes at render time. Do not
  change `bin/annotate_graph_svg.py`'s output: `diagnostics/graph.svg`
  must stay byte-identical (task 53 §3.1 relies on this). The fullscreen
  modal view must still scale up to fill the modal.

### 5.3 *Genome annotation* section

- Remove the "Genome map" `<h4>`. Move its info badge onto the
  *Genome annotation* `<h3>`.
- With the graph gone, the map, the "Annotation details" button and
  the GFF download take the full section width. Constrain the map to a
  sensible max width (for example, so a one-panel map is not wider than
  the viewport height) and centre it. The two-panel plastid map may use
  more width.
- Fix the legend clipping in `bin/organelle_map.py` (§2): compute the
  SVG height from the legend's real extent (three kind rows plus the
  footnote line plus a bottom margin), not from the fixed
  `LEGEND_ROW_HEIGHT * 4`. Add room for any new legend entry from §5.4.

### 5.4 Label barcode loci on the genome map

The map should name the barcode panel loci and **only** those, so they
stand out from the surrounding gene ring. The map already knows each
gene's name, so this is a name match against the panel config. No
barcode results are needed.

- **Input**: add a `--locus-panel` argument to `bin/organelle_map.py`,
  and pass it in `modules/local/organelle_map.nf` as
  `${file(params.locus_panel)}`, the same way
  `modules/local/extract_barcodes.nf` does (CONSTITUTION rule 13). Also
  pass `--assembly-target ${meta.assembly_target}` to select the panel.
  `main.nf` and the process's channel inputs do not change.
- **Matching**: a drawn CDS locus (after `cluster_loci`) gets a label
  when its primary `gene` is in the target's panel. Compare
  case-insensitively, exactly as `validate_barcodes.py` does
  (`.upper()` on both sides), so the map and `EXTRACT_BARCODES` cannot
  disagree about which genes are panel loci. Ignore `alt_genes`: an
  arc is labelled for what it was drawn as. tRNA and rRNA features are
  never labelled, because the panel is protein-coding only.
- **What a label means**: a label marks a *panel locus that was
  annotated*, not a *recovered barcode*. A panel gene that was
  annotated but then dropped by validation (e.g. below the identity
  floor) is still labelled. Pass, partial or dropped is the Barcodes
  tab's job (spec §6a.2 tab 4), and the map must not hint at it. The
  legend note must say "barcode panel loci", not "recovered barcodes".
  This is the honest reading of a name match (principle 7). The map
  has no access to validation outcomes, and showing them would mean
  wiring `EXTRACT_BARCODES` output into this stage. That is out of
  scope (§10).
- **Placement**: put the label radially *outside* the outer ring at
  the arc's mid-angle, with a short leader line. Anchor the text
  start/end by side of the circle so labels do not run back over the
  ring. Panels hold only a handful of barcode loci, so a simple
  collision nudge is enough (shift a label along the radius or angle
  when it would overlap the previous one). Do not build a general
  label-layout engine. Grow `PANEL_MARGIN` or compute it so labels fit
  inside the viewBox.
- **Plastid path2 panel**: labels follow the remapped features. A
  barcode feature dropped from path2 because it spans the SSC boundary
  gets no path2 label.
- Mark labelled features in the SVG (e.g. `data-barcode="1"` on the
  feature `<g>` and a `data-barcode-label` element). Tests can then
  assert on markup instead of geometry. Add a one-line legend note
  (e.g. "labelled: barcode panel loci").
- **Failure isolation is unchanged**: a panel file that is unreadable,
  or has no entry for the target, means no labels. It must never be an
  empty map or a non-zero exit. `main()` still never raises (the module
  docstring contract).

### 5.5 Bootstrap tooltips on the genome map

Task 53 gave the assembly graph Bootstrap-styled tooltips. The genome
map still uses the plain native SVG `<title>` tooltip, so the two
diagrams on the same tab behave differently. Apply task 53's helpers to
the map. This was recorded in `tasks/todo.md` (*Report UI backlog*) as
a task 53 follow-up, and is folded in here.

- **Markup is already compatible.** `_render_feature()` in
  `bin/organelle_map.py` emits each feature as `g[data-gene]` with its
  `<title>` as a **direct** child, which is the shape
  `promoteSvgTitleTooltips` requires. `organelle_map.py` needs no
  change for this, and the standalone `<id>.map.svg` keeps its native
  `<title>` tooltips, as task 53 §3.1 required for the graph.
- **Calls**: next to the existing
  `bindSvgModalMoveRestore('mapModal', 'organelle-map',
  'organelle-map-modal-body')`, add the same pair the graph uses:
  - `promoteSvgTitleTooltips('organelle-map', 'g[data-gene]')`
  - `hideSvgTooltipsAroundModal('organelle-map',
    'organelle-map-modal-body', 'g[data-gene]', 'mapModal')`

  The second call stops a tooltip being left floating when the
  fullscreen map modal opens or closes (task 53 §3.4).
- **Multi-line text.** Unlike the graph's one-line node text, the map's
  `_tooltip_text()` joins several lines with `\n`: gene,
  kind/coordinates/strand, source, `also called: …`, scores. A
  Bootstrap tooltip collapses newlines by default, so the text would run
  together. Keep the text plain (`html: false` stays, so a gene or
  contig name can never be read as markup). Instead, give the map's
  tooltips a custom class and set `white-space: pre-line` on that
  class's `.tooltip-inner` in `main.css`. Bootstrap 5.0.2 supports the
  `customClass` option, which can be set as a `data-bs-custom-class`
  attribute. The cleanest route is an optional `customClass` argument
  on `promoteSvgTitleTooltips` that sets that attribute. The graph call
  is unchanged and simply leaves it out. Check it works on 5.0.2 before
  relying on it. If it does not, scope the CSS by container instead.
  The tooltip `max-width` may also need widening for this class, since
  the coordinate line is long.
- **Barcode labels (§5.4)** are separate `text` elements with no
  `<title>`, so the helper leaves them alone. Hovering a labelled arc
  shows the arc's tooltip as usual.
- Delete the follow-up line from `tasks/todo.md` (*Report UI backlog*)
  once this lands, and the section heading too if it is then empty.

---

## 6. Barcodes tab

Render the *Panel barcodes not recovered* section only when
`barcodes_view.dropped` is non-empty. A `partial` is recovered, not
dropped (spec §6a.2 tab 4), so a sample with partials but no drops also
hides the section. The heading count "Recovered barcodes (n/N + k
partial)" already tells the reader that everything was found, so no
replacement line is needed. Apply the §4.4 genetic-code filter to this
tab's table.

---

## 7. Spec updates

Edit these in the same change. Follow the house rule: plain-text task
references, no links to task files, no line anchors.

- `spec/06a-reports.md` §6a.2:
  - Tab 1: Key findings shows severity badges (not row tint), sits
    beside the warnings mirror, and read QC uses a two-column layout.
  - Tab 2: gate decision is a headline badge, reads/bases/coverage are
    charts with exact values labelled, coverage floors are drawn on the
    chart, subsampling detail is in a modal, genetic codes are shown by
    name.
  - Tab 3: the graph sits under *Assembly statistics*, the assembly
    FASTA can be copied and downloaded, and the map labels barcode
    panel loci by name. State that a label means "panel locus,
    annotated", not "barcode recovered".
  - Tab 4: the not-recovered section shows only when non-empty.
- `spec/06a-reports.md` §6a.5 *Renderer inputs beyond `metadata.json`*:
  add the assembly FASTA.
- `spec/02-stages.md` stage 14 (`ORGANELLE_MAP`) entry: the stage now
  reads `params.locus_panel` to label panel loci, using the same name
  match as `EXTRACT_BARCODES`.

---

## 8. Work items

1. Overview: §3.1–3.4 (`overview.html`, `key-findings.html`, new badge
   macro, `read-qc.html`, `_top_blast_hit_finding` hint).
2. Validation: §4.1–4.4 (`validation.html`, `validation_view` if chart
   data needs shaping, genetic-code mapping and filter).
3. Assembly: §5.1–5.3 (`assembly.html`, `report.py` /
   `collate.py` / `render_report.py` FASTA plumbing, `main.css`,
   `organelle_map.py` legend height).
4. Genome-map labels: §5.4 (`organelle_map.py`,
   `modules/local/organelle_map.nf`).
5. Genome-map tooltips: §5.5 (`assembly.html`,
   `svg-title-tooltips.js` optional `customClass`, `main.css`, and
   removing the `tasks/todo.md` line).
6. Barcodes: §6.
7. Spec: §7.
8. Re-render and visually check (§9.3), then record the results in an
   Outcomes section on this task.

---

## 9. Tests

### 9.1 Unit — `scripts/tests/` via `scripts/pytest.sh`

All touched `bin/*.py` code needs 100% branch coverage (CONSTITUTION
rule 14). Then run flake8 on the touched `bin/*.py` and
`scripts/tests/*.py`.

`test_report.py`:

- Key findings: each of the five severities renders its badge (class
  and icon), the badge comes before the value text, and no
  `table-{class}` row tint remains. Update the existing `table-info`
  low-coverage assertion to the `info` badge. `low_coverage` renders no
  warning or danger badge on its Outcome row.
- The Top BLAST hit finding carries the compact hint, and other
  findings do not.
- Genetic-code filter: known ID, unknown ID ("NCBI table N"), `None`
  (`-`). It is applied in the Validation and Barcodes templates.
- Subsampling: modal and button present when `estimate.subsampled`,
  absent otherwise, and the "Subsampled: No" row present.
- Validation charts: chart containers render when their data is
  present, and a missing-input case renders "not available" rather
  than a chart with a zero bar. The coverage-floor values reach the
  chart script.
- Assembly FASTA: `assembly_fasta_src` is set for a real file and is
  `None` for absent or zero-byte files. The buttons render only when
  it is set. `render()` accepts the new argument.
- Graph placement: the `#assembly-graph` markup sits inside the
  *Assembly statistics* section, and the colour key sits outside
  `#assembly-graph`.
- Barcodes: the not-recovered section is absent when `dropped` is
  empty (including the partial-only case) and present when it is not.
- Map tooltips: when a map SVG is present, the rendered page calls
  `promoteSvgTitleTooltips` and `hideSvgTooltipsAroundModal` for
  `organelle-map` / `mapModal`. Neither is called when
  `organelle_map_svg` is absent. The server-rendered map markup still
  contains its `<title>` children, because conversion happens in the
  browser (same check task 53 §5.2 made for the graph).
- Existing self-contained / no-external-asset tests still pass.

`test_collate.py`: `target_fasta` is forwarded to the renderer.

`test_organelle_map.py`:

- Only CDS loci whose primary gene is in the target's panel get a
  label. A non-panel CDS gets none, and neither does a tRNA/rRNA.
- Case-insensitive match: a panel `rbcL` labels a drawn `RBCL`, and
  the reverse.
- The panel is selected per `assembly_target`: an `animal_mt` panel
  gene name does not label a `plant_pt` map.
- A clustered arc whose *primary* gene is in the panel gets exactly one
  label. An arc whose panel gene appears only in `alt_genes` gets none.
- Plastid path2: labels follow remapped features, and a feature that
  spans the SSC boundary has no path2 label.
- An unreadable panel file, or one with no entry for the target, gives
  a valid map with no labels and exit 0.
- **Legend containment**: every legend and label text `y` (plus its
  font size) lies within the viewBox height. This regression-tests the
  §2 clipping bug. Every label `x` lies within the viewBox width.

### 9.2 Integration — `tests/integration/assertions.sh`

- No new fixtures. The existing samples cover every branch: an
  `animal_mt` with barcodes (`INT-ANIMAL-01`), the plastid two-panel
  map (`INT-PLANT-01-pt`), and the `low_coverage` sample
  (`INT-PLANT-01-mt`).
- Extend the `ORGANELLE_MAP` block (task 44 assertions): for each
  assembling sample, the map SVG has at least one
  `data-barcode-label` element, and every labelled gene is in that
  target's panel in `assets/loci.json`. On `INT-PLANT-01-pt`, both
  panels carry labels. Assert on panel membership, not a fixed count,
  so the test does not drift with the reference data (CONSTITUTION
  rule 19).
- Confirm that the report text the existing report checks grep for
  (warn-floor alert, plastid alert) is unchanged.
- Run `-stub-run` end-to-end to confirm the `ORGANELLE_MAP` script
  change.
- A full integration run is needed because §5.4 changes a pipeline
  stage's command line.

### 9.3 Manual visual check (record in Outcomes)

Re-render `INT-ANIMAL-01`, `INT-PLANT-01-pt` and `INT-PLANT-01-mt`
from the integration outputs, and compare each against
`reference-material/task-54-report-int-animal-01.html` at desktop
width and at a narrow width (below `md`):

1. Every item in §3–§6 is present and laid out as described.
2. The graph fits its column, the zoom modal still fills the screen,
   and the colour key shows in both views.
3. The genome-map legend is fully visible, and barcode labels are
   legible and do not overlap on the plastid two-panel map.
4. Copy-to-clipboard puts a valid FASTA on the clipboard, and the
   download saves the same content.
5. Map tooltips: repeat task 53 §5.1 steps 1–5 against
   `#organelle-map` / `#mapModal`, on `INT-ANIMAL-01` (which has tRNA,
   rRNA and CDS features) in Chrome and Firefox:
   - Hovering shows a Bootstrap tooltip with **one line per field**
     and no native tooltip.
   - No tooltip is left behind when the map modal opens or closes.
   - The graph's tooltips are unchanged.
   - Opening `<id>.map.svg` directly still shows the native tooltip.
6. Save the report (`save-report.js`), reopen it, and repeat 2–5.

Give the rendered report paths to the user (CLAUDE.md).

---

## 10. Explicit non-goals

- Any change to sample status, gate logic, thresholds, severities or
  which findings and warnings exist.
- The run-level `run-report.html`. The new badge macro is written to be
  reusable, but adopting it there is a separate change.
- Tooltips on the barcode label text itself (§5.5): the arc carries
  the tooltip.
- Keyboard focus on map features, which task 53 §3.7 also left out for
  graph nodes.
- Labelling non-barcode genes on the map, or general label
  collision avoidance.
- Showing barcode validation outcomes (pass / partial / dropped) on the
  map. That would need `EXTRACT_BARCODES` output wired into
  `ORGANELLE_MAP`. If it is ever wanted, scope it as its own task.
- Changing `bin/annotate_graph_svg.py` or `diagnostics/graph.svg`.
- Upgrading Bootstrap (stays at 5.0.2) or adding an icon font.

---

## 11. Acceptance criteria

1. Overview: Key findings and Warnings sit side by side at 8/4, findings
   carry severity badges to the left of their values (no row tint),
   `low_coverage` reads as informational, and the top BLAST hit is in a
   smaller font.
2. Read QC: table and vertical-bar chart sit side by side, and the
   column reads "Clean reads".
3. Validation: the gate decision is a headline badge. The reads, bases
   and coverage charts show exact values without hovering, and the
   coverage chart draws all three floors. Subsampling detail is in a
   modal, reachable only when subsampling happened, and "Subsampled:
   No" is visible otherwise. No genetic code renders as a bare ID.
4. Assembly: the FASTA can be copied and downloaded. The graph sits in
   *Assembly statistics*, fits its column, and still zooms. The
   "Genome map" heading is gone, and its info badge is on *Genome
   annotation*.
5. Genome map: no clipped text. Annotated barcode panel loci are
   labelled outside the ring, and nothing else is. The legend says
   "barcode panel loci". A sample with no panel genes annotated still
   gets its map. Map features show Bootstrap-styled, line-preserving
   tooltips in the report and leave none behind around the modal. The
   standalone map SVG keeps native tooltips.
6. Barcodes: no "not recovered" section when nothing was dropped.
7. The report remains a single self-contained file. Unit tests pass at
   100% branch coverage on touched code, flake8 is clean, and stub-run
   and integration pass.
8. The spec is updated per §7.
