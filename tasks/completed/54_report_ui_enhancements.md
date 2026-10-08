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

---

## Outcomes (2026-10-07)

All eight work items in §8 were completed. Summary by area:

**§3 Overview** — `overview.html` puts Key findings (`col-md-8`) and
Warnings (`col-md-4`) in one row, Key findings first. Row tint on
`key-findings.html` was replaced with a new reusable
`macros/severity-badge.html` (`severity_badge(severity)`), rendered
before the value text on every row, carrying `title`/`aria-label` text
naming the severity in words. `_top_blast_hit_finding()` in
`bin/report/report.py` sets `'compact': True`; the template applies
`class="small"` from that flag only. `read-qc.html` is now a two-column
`col-md-6`/`col-md-6` layout, the "Clean reads (post CHOPPER +
FILTLONG)" header is now "Clean reads", and the filter-yield chart is
vertical (categories on x, retained % as a bar text label).

**§4 Validation** — the gate decision is now a bordered card with a
`fs-5` badge and a one-line meaning at the top of the section (the old
table row was removed). Three charts (`reads-chart`, `bases-chart`,
`coverage-chart`) replace the six reads/bases/coverage table rows, via
a new `renderNestedBarChart()` helper
(`scripts/report/static/js/nested-bar-chart.js`) and
`report.py::_validation_charts()` / `validation_view()['charts']`,
which carries `None`/empty values rather than a fabricated zero when
an input is absent (principle 7) — the JS falls back to a muted "Not
available" line when a chart has no traces at all. The coverage chart
draws the hard-min/warn/max floors as dashed lines + annotations, and
draws one narrow green bar per assembled target contig when
`flye_depth` has more than one. Subsampling detail moved into
`#subsamplingModal`, reachable from a button next to a
"Subsampled: Yes/No" row that always renders. Genetic codes render by
name via a new `macros/genetic-code.html` (`render_genetic_code`)
backed by `report.py`'s `GENETIC_CODE_NAMES` map and two new Jinja
filters (`genetic_code_name`, `genetic_code_tooltip`), applied in both
Validation and Barcodes.

**§5 Assembly** — `render()`/`build_context()` take a new
`assembly_fasta` argument, encoded as `assembly_fasta_src` via the
existing `_file_src()` helper (absent/zero-byte -> `None`, same
contract as `annotation_gff`). `bin/collate.py` passes
`args.target_fasta` as `assembly_fasta`; `bin/render_report.py` gained
a matching `--assembly-fasta` flag. *Assembly statistics* now shows a
copy-to-clipboard button (decodes the existing `data:` URI client-side
via `atob`, no second inlining) and a download link, both tooltipped
as "target contig(s) only, not the whole metaFlye assembly". The
assembly graph (`#assembly-graph`, its colour key, `#graphModal`,
`bindSvgModalMoveRestore`/`promoteSvgTitleTooltips`/
`hideSvgTooltipsAroundModal` calls) moved into a `col-sm-3` inside
*Assembly statistics*, alongside the stats table (`col-sm-4`) and
coverage chart (`col-sm-5`). The BandageNG SVG's overflow (fixed
`width="282.222mm" height="211.667mm"` root attributes, confirmed on
the real `INT-ANIMAL-01` graph) is fixed with
`#assembly-graph svg { max-width: 100%; height: auto; }` in
`main.css` — `bin/annotate_graph_svg.py` and `diagnostics/graph.svg`
are untouched (confirmed via `git diff --stat`). The "Genome map" `<h4>`
is gone; its info badge moved onto the *Genome annotation* `<h3>`
(now also describing what a barcode label means). The map's wrapper is
capped at `max-width: 1100px` and centred (`mx-auto`) — see Deviations
below.

**§5.4 genome-map labels** — `bin/organelle_map.py`: `load_locus_panel()`
loads and upper-cases `assets/loci.json`'s per-target panel (any
failure -> empty set, never raises); `cluster_loci()` takes `panel_genes`
and sets `is_barcode` on each locus from its *primary* gene only
(never `alt_genes`), and only for CDS loci (tRNA/rRNA from
`other_feats` are always `is_barcode: False`). `_render_feature()` adds
`data-barcode="1"`; a new `_render_barcode_labels()` draws each
labelled locus's gene name outside the ring at its mid-angle with a
leader line, anchored start/end by which half of the circle it's on, with
a simple radial nudge (`LABEL_NUDGE`) when two labels' angles are too
close. `render_svg()` grows the panel margin
(`LABEL_MARGIN_EXTRA = 90`) and the legend
(`_legend_extent()`/`_render_legend()` now compute extent from actual
content instead of a fixed `LEGEND_ROW_HEIGHT * 4`, fixing the §2
clipping bug for every map, not just labelled ones) when any locus is
labelled, and adds a "labelled: barcode panel loci" legend line.
`remap_to_path2()` carries `is_barcode` through unchanged (a
boundary-spanning feature is already dropped by the existing logic).
`--locus-panel`/`--assembly-target` were added to
`organelle_map.py`'s CLI and to `modules/local/organelle_map.nf`'s
script block (inline `${file(params.locus_panel)}`, same pattern as
`extract_barcodes.nf`); `main.nf`'s channel wiring is unchanged.
Verified directly against the real `INT-ANIMAL-01` (5 labels) and
`INT-PLANT-01-pt` (10 labels across both panels) integration outputs:
every legend/label text position lies inside the viewBox (script-
checked, see Testing below).

**§5.5 genome-map tooltips** — `assembly.html` adds
`promoteSvgTitleTooltips('organelle-map', 'g[data-gene]', 'map-tooltip')`
and the matching `hideSvgTooltipsAroundModal(...)` call next to the
existing `bindSvgModalMoveRestore('mapModal', ...)`.
`svg-title-tooltips.js`'s `promoteSvgTitleTooltips()` gained an
optional third `customClass` argument, setting
`data-bs-custom-class` (the graph's call omits it, unchanged). Verified
directly against the vendored `bootstrap.bundle.min.js`: it is
genuinely Bootstrap v5.0.2 (checked `Bootstrap v5.0.2` banner string)
and its `Tooltip` config schema already lists `customClass`, read
generically off any `data-bs-*` dataset attribute by `Manipulator.
getDataAttributes()` — so `data-bs-custom-class` → `customClass` works
on this exact pinned version without a library change. `main.css`
adds `.map-tooltip .tooltip-inner { white-space: pre-line; max-width:
320px; text-align: left; }`. `organelle_map.py` needed no change for
this (its `<title>` markup was already in the shape the helper
expects); the standalone `<id>.map.svg` is unaffected since the
conversion is client-side only in `report.html`.

**§6 Barcodes** — the "Panel barcodes not recovered" section in
`barcodes.html` is now wrapped in `{% if barcodes_view.dropped %}`;
the empty-row fallback text was removed since the heading count already
tells the reader everything recovered. The genetic-code filter was
also applied to this tab's table.

**§7 Spec** — `spec/06a-reports.md` §6a.2 (all four tab descriptions)
and §6a.5 (renderer-inputs bullet) updated; `spec/02-stages.md` stage
14 entry updated. All references are plain-text ("task
54_report_ui_enhancements.md"), no links to task files, per the house
rule.

### Deviations / judgment calls (read this section)

1. **Genome-map wrapper width (§5.3).** The task allows "the two-panel
   plastid map may use more width" than a one-panel map. Rather than
   compute a width from panel count (which risks turning into the kind
   of general sizing logic rule 19 warns against for a one-off
   presentation tweak), I used a single fixed `max-width: 1100px`
   wrapper for both cases. A two-panel map's SVG has a wider intrinsic
   `viewBox` (confirmed: `INT-PLANT-01-pt`'s is 1440×830 vs a one-panel
   map's 720×830), so it naturally renders closer to the 1100px cap
   than a one-panel map does, without any per-case logic — but both are
   still capped at the same absolute width. This is a conservative
   simplification, not a bug; flagging it explicitly since the task
   worded this as "may," not "must."
2. **Coverage-chart category labelling when multiple target contigs
   exist.** The task says to draw "one narrow green bar per target
   contig," without fully specifying how the single grey "whole" bar
   should coexist with several green ones in one chart. I chose to
   repeat the same grey "whole" value once per contig category (so it
   reads as a constant backdrop line across all green bars) rather
   than, e.g., drawing it as a single reference line spanning the
   chart. This path *is* exercised by this run's real data:
   `INT-PLANT-01-mt` (`low_coverage`, the `emit: all` multi-contig
   path) selected four target contigs
   (`contig_3`/`contig_5`/`contig_4`/`contig_6`, coverages
   `7.0/7.0/5.0/5.0`×), and its real `report.html` carries exactly the
   expected chart data: `{"labels": ["contig_3", "contig_5",
   "contig_4", "contig_6"], "values": [7.0, 7.0, 5.0, 5.0], "whole":
   28.48, ...}`. I did not view this rendered live in a browser (no
   browser tool available), so the choice of repeating the grey bar
   per category is confirmed correct in data but not confirmed to
   *read well visually* — worth a sighted look.
3. **"Subsampled" row literal wording.** The task's acceptance criteria
   quote `"Subsampled: No"` literally; I render exactly that string
   (and `"Subsampled: Yes"`) as one table-cell label, rather than
   splitting "Subsampled" into the row's label column and "No"/"Yes"
   into the value column (which the old row layout used for its other
   rows). This keeps the exact literal string intact and auditable
   without hovering, at a small cost to the row's visual consistency
   with its neighbours.

### Testing

- **Unit** (`scripts/pytest.sh scripts/tests/`): 748 passed, 100%
  statement+branch coverage across every file under `bin/`, including
  `bin/organelle_map.py`, `bin/report/report.py`, `bin/collate.py`.
- **flake8** (`/home/cameron/.local/envs/claude/bin/flake8`) clean on
  every touched `bin/*.py` and `scripts/tests/*.py`.
- **Stub-run**: `nextflow run . -profile stub,docker -stub-run`
  completed successfully (confirms `ORGANELLE_MAP`'s new
  `--locus-panel`/`--assembly-target` command line is well-formed).
  Needed `NXF_VER=24.10.0` — this sandbox's default `nextflow`
  (26.04.6) fails to even parse the pre-existing, untouched `main.nf`
  (a `error(...)` call's `+`-continued string literal), unrelated to
  this task; noting it here in case it affects other sessions'
  environment, not something this task fixes.
- **Integration**: `nextflow run . -profile integration,docker`
  (`NXF_VER=24.10.0`) completed, 66/66 tasks succeeded in 30m49s.
  `bash tests/integration/assertions.sh` then passed every assertion,
  including the new barcode-panel-label checks added to this file
  (every labelled gene on `INT-ANIMAL-01`'s and `INT-PLANT-01-pt`'s
  maps is a genuine member of that sample's own panel in
  `assets/loci.json`, and `INT-PLANT-01-pt`'s path1+path2 both carry
  labels).
- Additional structural checks run directly against the real
  integration output (not just fixtures) during this task, beyond
  `assertions.sh`: (a) decoded the `assembly_fasta_src` `data:` URI
  embedded in `INT-ANIMAL-01/report.html` and confirmed it is
  byte-identical to `organelle_assembly.fasta` on disk; (b) scanned
  every legend/label `<text>` position in all three samples' real
  `diagnostics/organelle_map.svg` files and confirmed none exceeds the
  viewBox (the §5.3 regression this task exists to fix); (c) confirmed
  `#assembly-graph svg` CSS is present in the rendered report and
  targets the real cause of the overflow (BandageNG's SVG root carries
  explicit `width="282.222mm" height="211.667mm"`).
- **Manual visual check (§9.3)** — rendered reports are at
  `tests/integration/output/INT-ANIMAL-01/report.html`,
  `tests/integration/output/INT-PLANT-01-pt/report.html`, and
  `tests/integration/output/INT-PLANT-01-mt/report.html`. No browser
  automation tool was available in this sandbox, so §9.3's checks were
  verified as follows:
  - Items 1 (layout per §3–§6), 2 (graph fits its column — CSS
    confirmed present and correctly targeted, colour key markup
    confirmed outside `#assembly-graph`), 3 (legend/label containment),
    4 (FASTA copy decodes/matches the real file) were verified
    **structurally** — either by direct inspection of the rendered
    HTML/SVG/JS or by extracting and checking embedded data
    programmatically (see above) — not by rendering in and interacting
    with a real browser.
  - Item 5 (live hover/tooltip behaviour, Bootstrap-vs-native, nothing
    left behind around the modal, in Chrome and Firefox) and item 6
    (save-report.js round trip) were **not** verified interactively —
    no browser tool was available. What *was* verified: the exact
    Bootstrap build vendored in this report (`v5.0.2`) supports the
    `customClass` tooltip option generically off `data-bs-custom-class`
    (confirmed by reading `bootstrap.bundle.min.js`'s own `Config`
    schema and `Manipulator.getDataAttributes()`, see Deviations/§5.5
    above), and that the JS calls wiring it up are present, correctly
    ordered relative to `bootstrap.Tooltip` construction in
    `index.html`, and covered by unit tests
    (`TestMapTooltips` in `scripts/tests/test_report.py`). A sighted
    reviewer with a browser should still complete item 5's five-step
    check and item 6 before treating this as fully verified — flagging
    this explicitly rather than claiming interactive verification that
    did not happen.


## Follow-up 2026-10-07

Looks OK but a few revisions required:

**Overview**
- Top BLAST Hit icon should not be <small>
- Input/run context row: first column should set the table to fill all available horizontal width, and second column should be col-auto, such that the first column+table consume all available horizontal width
- "Filter yield" chart - remove the heading, it's an unneccessary duplication of content. And make sure the bar chart fits the horzontal width. On my screen it's overflowing horizontally by about 120px. This should be a skinny bar chart - 50px width per bar, and the label should be rounded to an int (72.79% -> 73%).

**Validation**
- Stacked bar charts are terrible (reads aligned/recruited + bases raw/recruited). Each should be a single skinny (50px) stack composed of two bars:
    - Reads aligned: first bar shows recruited (green), second bar shows aligned but not recruited, full stack shows aligned (recruited is a subset of aligned)
    - Bases recruited: first bar shows recruited (green), second bar shows not-recruited, full stack shows raw reads (recruited is a subset of raw reads)
- Estimated coverage chart is ok, but ridiculously wide and overflowing the column by about 600px. It also doesn't seem to be stacked - again, there should be two stacked bars, with full stack showing recruited bases, first bar showing assembled bases and second bar showing off-target bases. The threshold labels font should be ~2pt larger.
- Recruitment/coverage gate table is redundant - floors are shown in bar chart, and "Subsampled" doesn't need to be in a table. Suggest moving to a standalone "Subsampling" button which is disabled with tooltip "No subsampling done" when no subsampling was done.

**Assembly**
- Looks good, but Download button should use a download icon, not text, and these buttons should stack horizontally with 1rem of margin between them.

### Follow-up resolution (2026-10-07)

All six points addressed:

- **Top BLAST hit**: `<small>` now wraps only `finding.text`, not the
  severity badge — the badge stays full-size (`key-findings.html`,
  `scripts/tests/test_report.py`'s
  `test_compact_hint_wraps_text_in_small_not_badge`).
- **Inputs / run context row**: `overview.html`'s first column is
  `col` (fills remaining width), second is `col-auto`; `inputs.html`'s
  table gained `w-100` so it actually stretches to fill that column
  rather than shrink-wrapping its content.
- **Filter yield**: the "Filter yield: X% ... Y% ..." sentence (which
  restated numbers the chart already showed as bar labels) is replaced
  with a plain caption that keeps only the CHOPPER + FILTLONG tool
  mention. Bar labels round to the nearest integer. Bars render at a
  computed width targeting ~50px regardless of column width.
- **Validation reads/bases/coverage charts**: `nested-bar-chart.js` was
  rewritten from Plotly `barmode: 'overlay'` (two separate bars, one
  narrower) to a genuine `barmode: 'stack'`: a green "subset" segment
  stacked under a grey "remainder" segment (`whole - subset`, clamped
  to 0), so the full stack height reads as the whole and the green
  segment as its share — matching the brief's "first bar
  recruited/assembled, second bar the remainder, full stack the
  whole" description exactly, for reads, bases and coverage alike.
  Floor-line annotation font raised 9 → 11. Bars target ~50px,
  computed from the container's actual rendered width.
- **Chart overflow (120px / 600px)**: root cause was Plotly rendering
  at a fallback width before the container had a real size — on
  Overview this is a transient layout-timing issue, and on Validation
  it is the hidden `tab-pane` (`display: none` until its tab is
  clicked) reporting `clientWidth` 0 at render time. Fixed uniformly
  via a new `renderWhenSized()` helper in `nested-bar-chart.js`: if the
  container has no width yet, it defers the `Plotly.newPlot` call
  (and the pixel-width bar-sizing calculation) to the first
  `ResizeObserver` callback that reports a non-zero width, then keeps
  observing via `Plotly.Plots.resize()` for later container changes —
  the same pattern already used for the assembly coverage chart
  (`assembly.html`, predating this task). The filter-yield chart on
  Overview now also calls `renderNestedBarChart()` instead of its own
  bespoke Plotly script, so it gets the same treatment.
- **Gate table / Subsampled row**: the "Fail / warn / max coverage
  floors" row is removed (redundant with the coverage chart's floor
  lines). The sibling-panel rows (target-assigned bases,
  sibling-assigned bases, sibling fraction) are now their own
  conditional table, shown only when that data exists — previously
  they shared a table with the two removed rows, so an empty table
  would otherwise render when sibling data is absent. "Subsampled" is
  now a standalone button below that table: a clickable "Subsampling"
  button (opens the existing modal) when
  `validation_view.estimate.subsampled`, otherwise a disabled
  "Subsampling: none" button inside a Bootstrap-tooltip wrapper
  reading "No subsampling done" (`data-bs-toggle="tooltip"` on a
  wrapping `<span>`, since disabled `<button>` elements don't fire the
  hover events Bootstrap's tooltip needs). The button's own text
  ("Subsampling" vs "Subsampling: none") still states the fact
  without requiring a hover, for CONSTITUTION rule 18 (a printed
  report has no hover) — the tooltip is a bonus for the interactive
  case, not the only way to read the state.
  `scripts/tests/test_report.py`'s `TestSubsamplingModal` tests were
  updated to match (the literal "Subsampled: Yes"/"Subsampled: No"
  strings this task's original work introduced are superseded by this
  follow-up).
- **Assembly FASTA buttons**: the copy button already had an inline
  SVG icon (unchanged); the download `<a>` now carries a Font Awesome
  "download" SVG glyph (512×512 viewBox, licensed the same as the
  file's existing icons) instead of the word "Download", and both
  controls sit in a `d-flex` wrapper with `gap: 1rem` so they lay out
  horizontally with the requested spacing.

Verification: re-rendered `INT-ANIMAL-01` (subsampled, single
contig), `INT-PLANT-01-mt` (not subsampled, 4-contig multi-target) and
`INT-PLANT-01-pt` via `bin/render_report.py` inside the
`neoformit/daff-wf5-scripts:test` image against the existing
integration output fixtures, and inspected the rendered HTML/JS
directly for each change above (badge/`<small>` nesting, `w-100` +
`col-auto`, the filter-yield caption and rounding, the stacked-chart
JS and its `remainderName`/`subsetName` options including the
multi-contig coverage case, both subsampling-button states, and the
download icon + `gap: 1rem` wrapper). `scripts/pytest.sh scripts/tests/`
passes all 748 tests at 100% branch coverage on every touched
`bin/*.py` file, and flake8 is clean on all touched `bin/*.py` and
`scripts/tests/*.py`. No interactive browser check was done for this
follow-up round either (same constraint as the original task's
Outcomes note) — the resize/tooltip/stack behaviour is verified via
the rendered markup and JS, not a live hover/click pass.


## Follow-up 2026-10-08

### Validation
**Bar charts**
- Please remove value labels. Value on hover is enough.
- Hover labels are far too busy. Is it possibly to show just the value only?
- Estimated coverage chart should not be stacked (I was wrong about this one,
  sorry!). Please change to two standard side-by-side bars. The floor/ceil
  label font is still too small, please bump another 2px.

**Annotation cross-checks**
- Table should be small text please
- Table badges:
  - Agreed (both annotators): Shows "Raw", "Canonical" - is this a bug? If not,
    what does this mean?
  - Miniprot only: as above
  - Non-CDS annotator only: "lagli (off-panel)" - should have brief tooltip to
    explain what this mean. Is lagli a gene?

**New table**
This section should start with a table of two columns and three rows that will absorb some components that are currently floating around the section:
- Status: this should include the "OK coverage cleared both the fail and warn
  floors"
- Subsampling: if yes, show "subsampling" button renamed to "Details"; if no,
  say "Below threshold - not performed"
- Warnings: if there are warnings, show them here. If not, "None".

### Assembly
- Looking good, but the "Assembly statistics" table should be small font

**Per-contig mean coverage bar chart**
- Aggregate off-target contigs into a single
  bar. Change display to vertical bars. Remove the legend - now replaced by the
  x-axis labels.
- Remove plotly title and replace with HTML title matching "Assembly graph".
  Move table caption into info badge consistent with "Assembly graph".

### Barcodes section
- Looking good, but we should include a "Download" icon button on each barcode
  row, and a "Download FASTA" button before the table, to download all barcode
  sequences.

### Follow-up resolution (2026-10-08)

**Bandage graph tooltip** (separate request, same day, addressed first):
trimmed the info-badge tooltip on "Assembly graph" to 2-3 sentences; the
full explanation (hover semantics, colour-key meaning, mixed/unknown
node handling) moved to a caption below the colour key inside the
fullscreen graph modal (`assembly.html`).

While verifying that, found `scripts/tests/test_report.py`'s
`TestAnnotationDetailsModal` asserting on the stale button text
"Annotation details" against a template that (outside this session)
already reads "Annotation statistics" — updated the test to match,
not the template.

**Barcodes**: "Download barcodes FASTA" moved above the table. Each
recovered-locus row gets an icon-only download button next to the
existing "View" button (`macros/sequence-display.html`'s new
`render_download_sequence_button`, `downloadSequenceAsFasta()` added
to `sequence-display-modal.html`'s shared script) — builds the FASTA
client-side from the sequence already inlined for the View button, so
nothing new is embedded in the page.

**`organelle_map.py` legend note removed**: per explicit instruction,
not patched — `include_barcode_note` was deleted outright from
`_legend_extent()`/`_render_legend()`/`render_svg()` (it had been left
mid-edit outside this session, emitting an unclosed `<text>` tag with
no content, which broke 2 tests). The feature this supported — a
"labelled: barcode panel loci" line in the map legend — is gone; the
barcode labels themselves (§5.4, `data-barcode`/`data-barcode-label`
attributes, the leader lines) are untouched and still work, since they
never depended on this legend note. 3 tests covering the removed
behaviour were deleted/trimmed from `test_organelle_map.py`.

**Validation tab — bar charts**: removed always-on bar text labels
from `nested-bar-chart.js` across all three charts (reads, bases,
coverage) and the Overview filter-yield chart that shares the same
helper — hover now shows just the bare value (`hovertext` +
`hoverinfo: 'text'`, no trace name/category clutter). This is a
deliberate, explicit exception to CONSTITUTION rule 18's usual
"visible without hovering" default for these particular charts, per
direct instruction — flagging it here since it departs from this
task's own stated principle rather than silently following it. The
coverage chart changed from stacked to `barmode: 'group'` (two
ordinary side-by-side bars: "Recruited pool estimate" and "Assembled
(Flye)") — added an `opts.mode` ('stack' default, 'group' opt-in) to
`renderNestedBarChart()` rather than a second function, since the
pixel-width/resize/hover machinery is identical either way. Floor-line
label font raised 11 → 13 (this task's own first follow-up round had
already raised it 9 → 11).

**Validation tab — Annotation cross-checks**: fixed a genuine bug, not
just a style request. `cds_crosscheck.agreed`/`miniprot_only` can be a
flat raw-name list (`annotate_summary.cds_crosscheck()`'s own output)
*or* a `{'raw': [...], 'canonical': [...]}` pair once `collate.py`'s
`with_canonical_names()` has paired them up for real pipeline runs
(§5.3) — the template was iterating `crosscheck.agreed` directly,
which for the dict shape iterates its *keys* ("raw", "canonical"),
rendering those literal strings as badges instead of the gene names
inside them. Fixed in Python (`report.py`'s new `_crosscheck_view()` /
`_crosscheck_names()` / `_crosscheck_annotator_only()`, keeping
severity/shape logic out of the template per this task's own rule),
which normalises both shapes into `{'display', 'raw'}` dicts — the
canonical name badge now carries an "also called: <raw>" tooltip only
when the two differ. `annotator_only` entries (`{'gene', 'reason'}`,
optionally `+ 'canonical'`) get the same `display`/`raw` treatment
plus a `reason_text` tooltip from a new `ANNOTATOR_ONLY_REASONS` dict
(`off_panel` / `overlap` / `no_exon_data`) — answering the question
asked directly: yes, `lagli` is a real gene (a mitochondrial
intron-encoded LAGLIDADG homing endonuclease), and `off_panel` means
MITOS2 called it but it isn't one of the target's protein-coding
barcode/annotation genes, so miniprot was never going to call it
regardless of coordinates (this is exactly the case
`rescue_disqualification()`'s own docstring documents for
`INT-ANIMAL-01`). Table rendered `small`.

**Validation tab — new Status/Subsampling/Warnings table**: replaced
the floating headline-badge card, the two hardcoded alert blocks
(low_coverage / sibling-split-unavailable), and the standalone
Subsampling button with one `table-sm` table of three rows, absorbing
all of it: **Status** (badge + meaning text + info badge, unchanged
content), **Subsampling** (the existing modal-opening button, relabelled
"Details", when subsampled; literal "Below threshold — not performed"
otherwise — the disabled-button-with-tooltip approach from this task's
first follow-up round is gone, since the request here was a plain
text fallback instead), **Warnings** (both alert texts verbatim when
either condition holds, "None" otherwise). Keeping the exact existing
wording for the low_coverage alert (rather than switching to the
similar-but-differently-worded text in the global `warnings` list used
by the Overview mirror) was a deliberate choice: `tests/integration/
assertions.sh` greps `INT-PLANT-01-mt/report.html` for "below the
coverage warn floor" and `test_report.py` pins the same substring
appearing in both the Overview and Validation tab ranges — reusing the
global list's wording would have still satisfied both (the substring
is common to both versions), but keeping the Validation tab's own
literal text avoids relying on that coincidence.

**Assembly tab**: "Assembly statistics" table rendered `small`.

**Verification**: `scripts/pytest.sh scripts/tests/` — 750 tests pass,
100% branch coverage on every touched `bin/*.py` file (`report.py`,
`organelle_map.py`). flake8 clean on all touched `bin/*.py` and
`scripts/tests/*.py`. Re-rendered `INT-ANIMAL-01` directly via
`report.render()` (passing `barcodes_fasta` explicitly, since
`bin/render_report.py`'s CLI doesn't expose that flag — a pre-existing
gap, not touched here) against the real integration fixture and
inspected the output HTML/JS for every change above: the Status table,
the "Details"/"Below threshold" subsampling states, the corrected
ATP6/ATP8/.../LAGLI badges with their tooltips (confirming the bug fix
against real `INT-ANIMAL-01` data, not just a synthetic test), the
`mode: 'group'` coverage chart call, the absence of bar-text-label code
in `nested-bar-chart.js` specifically (a `textposition` string
elsewhere in the page is the unrelated, pre-existing per-contig
coverage chart in `assembly.html`), font-size 13 floor labels, and the
`small` classes on both tables. No interactive browser check was done
(same constraint noted in every previous round of this task) — hover
behaviour and visual bar width/spacing in `group` mode are verified via
the rendered JS config, not a live hover/screenshot pass.
Updated in place: `tests/integration/output/INT-ANIMAL-01/report.html`.


## Follow-up 2026-10-08 (2)

### Validation chart data audit

#### Findings

A review of the Validation tab's "Read recruitment" charts traced each
bar back to its pipeline source and found they no longer meant what
their labels said:

| Chart | Bar | Source | Actual meaning |
|---|---|---|---|
| Reads | Recruited | `recruit_stats.json` `reads_recruited` | Reads passing `recruit_filter.py` thresholds |
| Reads | Not recruited | `reads_aligned − reads_recruited` | Reads minimap2 aligned to the target panel that the filter then dropped |
| Bases | Recruited | `sample_status.json` `total_recruited_bases` | All bases in the recruited FASTQ |
| Bases | Not recruited | NanoPlot **raw** bases − recruited bases | QC losses and non-recruitment mixed together |

Pipeline order: RAW → CHOPPER → FILTLONG (clean) → RECRUIT (align to
target panel, then filter) → COVERAGE_GATE.

- The Reads chart carried no information: `recruit_thresholds` ships
  at 0/0 (task 28 §10), so `reads_recruited == reads_aligned` on every
  fixture and "Not recruited" was always 0. Its denominator was also
  aligned reads, not the whole sample.
- The Bases chart spanned two pipeline stages. On `INT-ANIMAL-01`,
  raw 7.90 Mb → clean 5.76 Mb → recruited 5.57 Mb, so ~92% of the
  "Not recruited" bar was really QC loss.
- The two charts used different denominators, so they could not be
  read together.

The root cause goes back to §4.2: it drew separate table rows ("Total
raw bases", "Total recruited bases") as whole vs subset, which implied
a subset relationship that crosses the QC stage.

#### Requested changes

- [x] Replace the Reads and Bases charts with a single bases chart with
  three bars: **Raw**, **Clean**, **Recruited**, in Mb. Each step is a
  true subset of the one before, and each drop maps to a real pipeline
  stage (QC, then recruitment).
- [x] Estimated coverage chart: all bars green. Add a **Subsampled** bar
  between Recruited and Assembled when subsampling was done.
- [x] All Validation bar charts: hover shows the value only, since the
  x-tick labels now name each bar.
- [x] Estimated coverage chart: log-scale y-axis (requested mid-change).

#### Resolution

- `report.py`: `_validation_charts()` now returns labelled bars only.
  `_bases_bars()` gives Raw / Clean / Recruited in Mb (`BASES_PER_MB`).
  `_coverage_bars()` gives Recruited (`estimated_cov`), Subsampled
  (`post_subsample_cov`, only when `estimate.subsampled`), then Flye
  depth: "Assembled" for one target contig, "Assembled (contig_N)" for
  several, so multi-contig samples keep their per-contig breakdown
  (§4.2). `post_subsample_cov` is `estimated_cov` × the realised
  subsample fraction (`coverage_gate.py`), so all three bars share a
  basis. The unused `recruitment` / `raw_bases` keys were dropped from
  `validation_view`.
- `nested-bar-chart.js`: new `renderLabelledBarChart()` (one bar per
  x tick, single colour, value-only hover, optional floors, optional
  `logY`) replaces `renderPairBarChart()`. `renderNestedBarChart()` is
  now only used by the Overview filter-yield chart, so its group mode,
  floors and pixel-target sizing were removed as dead code.
- `validation.html`: one "Bases (Mb)" chart and one "Estimated
  <organelle> coverage" chart side by side; info-badge text rewritten
  to describe what each bar actually measures.
- `spec/06a-reports.md` §6a.2 tab 2 updated to match.

Judgment calls:

- **Bases chart colour**: all green, matching the coverage chart. Only
  the coverage chart's colour was specified.
- **Log-scale y-axis is visible** (ticks and axis line) on the coverage
  chart, unlike the hidden linear y-axis elsewhere: on a log scale, bar
  heights are not proportional to the values, and a hidden axis would
  invite misreading them.
- **Overview filter-yield chart hover unchanged** ("Retained: 73%"):
  its stacked segments are not named by x ticks, so the segment name is
  still the only on-hover cue to which colour is which.

Verification: 752 tests pass, 100% branch coverage on every `bin/*.py`
file, flake8 clean, `nested-bar-chart.js` passes `node --check`.
Re-rendered `INT-ANIMAL-01` (subsampled, single contig: Recruited
327.52× / Subsampled 294.17× / Assembled 123.0×; bases 7.90 / 5.76 /
5.57 Mb), `INT-PLANT-01-mt` (four target contigs, not subsampled) and
`INT-PLANT-01-pt`, and checked the chart data in each. Not checked in
a live browser.

#### Follow-up: title and y-axis on the bases chart

- [x] Chart title "Bases (Mb)" -> "Read recruitment".
- [x] Enable the bases chart's y-axis (previously hidden on all
  linear-scale `renderLabelledBarChart` charts) with label
  "Megabases".

`nested-bar-chart.js`'s `renderLabelledBarChart()` gained
`opts.yAxisTitle`: when set on a linear (non-log) chart, it shows a
titled, `automargin`, solid-axis y-axis instead of the default hidden
one — other linear charts that don't pass it (none currently) keep
the old hidden-axis behaviour. `validation.html`'s `bases-chart` call
now passes `yAxisTitle: 'Megabases'` and its heading reads "Read
recruitment". Verified: `node --check` passes, 181/752 tests still
pass at 100% branch coverage, flake8 clean, and all three fixtures
re-rendered with the new title and `yAxisTitle: 'Megabases'` confirmed
in the embedded chart script.

**Bug found after the above**: the title didn't actually render — user
reported "no y-axis title shown". Root cause: in the pinned Plotly
3.0.0 bundle, an axis `title` is a nested object
(`title: {text, font, standoff}`, confirmed against the bundled
schema in `plotly-basic-3.0.0.min.js`), not a plain string; passing a
bare string is silently dropped by schema coercion rather than
raising. `title: opts.yAxisTitle` was wrong in three places —
`renderNestedBarChart`'s `yaxis.title` (used by the Overview
filter-yield chart's `%` label, also broken, not previously reported),
`renderLabelledBarChart`'s `yaxis.title` (new in this follow-up), and
its unused `xaxis.title` (no template passes `xAxisTitle` yet, but
same bug). All three fixed to `title: { text: ... }`. Re-verified:
`node --check` passes, all three fixtures re-rendered with
`title: { text: opts.yAxisTitle }` confirmed in the embedded script;
full suite still 752 tests / 100% branch coverage (no Python touched
by this fix). Not checked in a live browser.

### Genome annotation

Bandage adds this into the SVG which creates an annoying hover-effect which
conflicts with our nice Bootstrap tooltips:

`<title>Qt SVG Document</title>`

Please remove this from the SVG when rendering into the report.

- [x] Strip the root `<title>Qt SVG Document</title>` from the
  BandageNG SVG before it's embedded in the report.

#### Resolution

`bin/annotate_graph_svg.py`'s `strip_xml_prolog()` (already run on
every SVG in `run_annotate()` before the per-node `<title>` tooltips
are added) now also strips that exact root title via a new
`_ROOT_TITLE_RE`. It only matches the literal Qt boilerplate text, so
the per-node tooltip `<title>` elements this module adds later are
untouched. Verified: added `test_root_qt_title_stripped` and
`test_per_node_title_left_alone` to `StripXmlPrologTests`; 51/51 tests
pass in `test_annotate_graph_svg.py` with 100% branch coverage on
`annotate_graph_svg.py`; full suite (752+ tests) still 100% branch
coverage throughout; flake8 clean. Not re-rendered against a real
BandageNG SVG fixture (none of the checked-in fixtures carry the real
Qt title text), so the fix is verified at the regex/unit level only.


## Follow-up 2026-10-08 (3)

### Overview

**Run context table**
- Shows nothing for started, completed, walltime. I would expect this to be
  displayed for an integration test run?
- Move "Bundle" table into the "View tool versions" modal, under a subheading
  "Refernce data bundle"

**Input param buttons**
- Rename "View all parameters" to "Run parameters"
- Rename "View tool versions" to "Tool versions"

### Validation

**"Estimated coverage" chart**
- The log-scale Y-axis is strange, as it shows "2, 5" as minor tick labels
  between "10, 100". Suggest removing minor tick labels.
- If possible, please add a hover label to the floor/ceil labels showing the
  value.

- [x] Run context table blank fields.
- [x] Move "Bundle" table into the "View tool versions" modal, under
  "Reference data bundle".
- [x] Rename "View all parameters" -> "Run parameters", "View tool
  versions" -> "Tool versions".
- [x] Estimated coverage chart: remove log-axis minor ticks.
- [x] Estimated coverage chart: hover on floor/ceiling lines shows the
  value.

#### Resolution

- **Run context table blanks**: not a pipeline bug. `wall_time`
  (`report.py`'s `wall_time_context()`) is populated from a
  `workflow_start` argument that only `render()`'s caller supplies —
  `collate.nf`/`collate.py` pass `workflow.start` from Nextflow on a
  real run. The blank fields the user saw came from my own earlier ad
  hoc `report.render()` calls in this session (used to re-check the
  chart work), which called `render()` directly and left
  `workflow_start` unset. Re-rendered all three fixtures with an
  explicit `workflow_start` to confirm the field populates correctly
  — no code change.
- `inputs.html`: buttons renamed to "Run parameters" / "Tool
  versions" (modal titles unchanged — not requested). The "Reference
  data bundle" table (Bundle, Pipeline commit, Reference bundle,
  Bundle built) moved from `provenance.html` into the tool-versions
  modal body, under a `<h5>Reference data bundle</h5>` heading.
  `provenance.html` is now unused (its only content moved) and was
  deleted; `overview.html` no longer includes it.
- `nested-bar-chart.js`: log-scale y-axis now sets `dtick: 1`, which
  restricts ticks to whole powers of ten and removes the unlabelled-
  looking 2/5 minor ticks Plotly adds by default between decades.
  Floor/ceiling lines were a layout `shape` (not hoverable in Plotly)
  with a separate always-visible label annotation; they're now a
  `scatter`/`mode: 'lines'` trace per floor (still dashed, same
  colour) with `hovertext` set to `"<label>: <value>"`, so hovering
  anywhere along the line shows the floor's name and value. The
  static annotation label at the chart's right edge is unchanged.

Verification: 754 tests pass (+2 vs. previous count — unrelated to
this follow-up, pre-existing), 100% branch coverage on every touched
`bin/*.py` file (no Python touched here), flake8 clean, `node --check`
clean. All three fixtures re-rendered; confirmed in the embedded
script: `Run parameters`/`Tool versions` button text, the "Reference
data bundle" table under the versions modal, populated
started/completed/walltime rows, `dtick: 1` on the log axis, and the
floor-line hover trace with label+value text. Not checked in a live
browser.

Two follow-up corrections after this:

- The `wall_time` fixture re-render above used a hardcoded
  `workflow_start` ("2026-10-08T09:00:00") that didn't account for the
  Docker container's clock being ~8 hours behind that guess, so
  "completed" came out before "started" and `wall_time_context()`'s
  `max(delta, timedelta(0))` clamp silently hid the negative duration
  as "00:00:00". Not a code bug — re-rendered all three fixtures with
  `workflow_start` computed from the container's own `datetime.now()`
  instead of a guessed literal.
- The floor-line hover added above "aren't really working" (user
  report) — reverted: floors are layout `shapes` again (not a
  hoverable trace), and `floor.label` now carries the value itself
  (e.g. "Fail 10x", "Warn 50x", "Max 300x") built in `validation.html`
  as `'Fail ' + coverage.floors.hard_min + 'x'` etc., since the line
  itself has no hover. `nested-bar-chart.js`'s `renderLabelledBarChart`
  doc comment updated to say `floor.label` should include the value.
  Verified: 181/754 tests pass, 100% coverage, flake8/`node --check`
  clean, all three fixtures re-rendered with the new label text
  confirmed in the embedded script (source only shows `'Fail '` etc.
  since the value is concatenated in-browser at render time).

- [x] Bump x-tick label font size on the bases/coverage charts.

`renderLabelledBarChart`'s shared `xaxis.tickfont` set to `{ size: 14 }`
(both charts use this one function). Verified: `node --check` clean,
754 tests pass, all three fixtures re-rendered.
