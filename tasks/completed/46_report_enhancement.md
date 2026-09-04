# Task 46 — Per-sample `report.html`: presentation review follow-ups

**Phase:** P4a (spec 06-phases.md).
**Goal:** Act on the project owner's review of task 43b's rendered
`INT-ANIMAL-01` report — two real defects, one tab reorder, and a set
of layout/presentation refinements across all four tabs.

**Depends on:** task 43b_report_stage_tabs.md, which shipped the
four-tab report this reviews.

**Related tasks:** 44_organelle_map.md supplies the real SVG for the
Assembly tab and is unaffected by this task. 45_run_report.md reuses
the shared components — **this task should land before it** (§13).

---

## 1. Overview

Task 43b shipped the four reader-perspective tabs. The project owner
reviewed the rendered `INT-ANIMAL-01` report and returned a list of
corrections. Most are presentation: reordering blocks, collapsing
adjacent tables, moving a dense section behind a modal.

Two are not presentation. Investigating the review surfaced a **badge
rendering defect that makes five elements invisible in every report**,
and established that **the report under review was rendered from stale
fixture data**, which explains two "blank section" observations that
are not code defects at all. Both are recorded in §2 before the
cosmetic work, because one of them changes what the review means.

Two review items conflicted with CONSTITUTION principle 7 and were
resolved with the project owner before this brief was written (§3).

---

## 2. Two findings that are not presentation

### 2.1 `text-bg-*` does not exist in the vendored Bootstrap — five badges render invisible

`scripts/report/static/css/bootstrap.min.css` is **Bootstrap v5.0.2**.
The `text-bg-*` contextual-badge helpers were introduced in **Bootstrap
5.2**. They are absent from the vendored stylesheet — `grep -c
'text-bg-success' bootstrap.min.css` returns 0.

`.badge` in 5.0.2 sets `color:#fff` and **no** background. So every
`text-bg-*` badge renders as white text on a transparent background:
invisible against the report's white page.

This is the actual cause of the review's *"Input params — status row is
blank"*. The status row is not empty; its badge is white-on-white.

Five call sites are affected, across three templates:

| Template | Badge |
|---|---|
| `components/inputs.html` | sample status (task 43a) |
| `components/validation.html` | gate decision |
| `components/validation.html` | "below species threshold" |
| `components/validation.html` | `annotator_only` cross-check genes |
| `components/barcodes.html` | barcode drop-out reason |

The review's own suggestion for the Barcodes tab — `class="badge
text-dark bg-warning"` — is exactly the right 5.0.2-compatible form.
Apply that pattern to **all five**, not just the one that was noticed.

**Decision: fix the call sites, do not upgrade Bootstrap.** Upgrading
the vendored stylesheet to 5.3 would touch every layout primitive in a
report nobody has visually regression-tested, to gain one helper class
that a two-line macro replaces. Rule 19 — the vendored asset is pinned
deliberately.

**Build a single Jinja macro** (e.g. `macros/badge.html`'s
`render_badge(text, context, title='')`) that owns the light/dark
pairing, rather than hand-pairing `bg-*` with `text-*` at each call
site. The pairing is not arbitrary — it follows the background's
luminance, and getting it wrong reintroduces an unreadable badge:

| Context | Background | Foreground |
|---|---|---|
| `success` | `bg-success` (#198754) | `text-white` |
| `danger` | `bg-danger` (#dc3545) | `text-white` |
| `secondary` | `bg-secondary` (#6c757d) | `text-white` |
| `warning` | `bg-warning` (#ffc107) | `text-dark` |
| `info` | `bg-info` (#0dcaf0) | `text-dark` |

One macro also means task 45_run_report.md inherits the fix instead of
forking the broken pattern (§13).

### 2.2 "Sequencing data quality is blank" is stale fixture data, not a defect

The review asks whether the empty read-QC section is an artefact of the
integration run. It is — but not in the way the question assumes.

`tests/integration/output/` is **gitignored**. The
`INT-ANIMAL-01/metadata.json` the reviewed report was rendered from has
an **mtime of 2026-09-02 12:54**, and carries **no `read_qc` key and no
`coverage.recruitment` key at all**. Both were added by task 43a, in
commit `1b5ec05`. The file predates the plumbing that populates them.

So the reviewed report was rendered from a pre-43a pipeline run. The
renderer behaved correctly: `read_qc: null` renders "No read QC data is
available for this sample", which is the honest output for that input
(rule 18 — a missing input must not be rendered as a zero result).

Two consequences the review could not have seen:

- The Validation tab's **"Reads aligned / recruited" row is silently
  absent** for the same reason — `coverage.recruitment` is missing, so
  the `{% if %}` guarding that row is false. Same root cause, same
  non-fix.
- **Nothing in this task should "fix" either section.** The action is
  to **re-run `-profile integration`** and re-review against fresh
  output. 43a's own outcomes record a real run producing populated
  `read_qc` (72.79 % reads / 72.92 % bases retained), so the plumbing
  is known-good.

**Do, however, verify the NanoPlot directory nesting during that
re-run.** The published tree shows
`qc/nanoplot_raw/nanoplot_raw/NanoStats.txt` — a doubled directory.
That doubling is a `publishDir` artefact (the module publishes to
`.../qc/nanoplot_raw` and emits a directory itself named
`nanoplot_raw`), and the path *staged into COLLATE* should be the inner
directory, so `parse_nanostats` should resolve `<dir>/NanoStats.txt`
correctly. Confirm that against a real run rather than assuming it;
if `read_qc` is still null on fresh output, the nesting is the first
place to look.

---

## 3. Two review items that conflicted with principle 7

Both were put to the project owner before this brief was written; both
were resolved in favour of the principle. Recorded here so the
resolution is not silently re-litigated during execution.

### 3.1 Key findings stays above the fold

The review asked for Input params + Run context at the top of Overview,
"followed by warnings and then key findings". Key findings last
contradicts spec §6a.1 ("Overview is tab 1 and opens by default, so the
first thing a reader sees is *Key findings*"), spec §6a.2 ("it is the
first thing on the first tab"), brief.md §3.8's "negative clarity
front-and-centre", and task 43b's
`test_key_findings_is_inside_overview_tab_not_a_header_row`.

**Resolved:** Input params and Run context become a **compact
two-column strip at the very top** of the Overview tab; **Key findings
sits immediately below it**; warnings below that; read QC last. The
review's intent — orient the reader before the findings — is served,
and negative clarity does not move below the fold. **No spec change
is required**, and 43b's DOM-order test still holds (§11).

Both panels are short (8 rows and ~9 rows), so a slim strip is
achievable without a scroll. Keep them visually subordinate — smaller
type, tight table styling — so Key findings remains the first thing
the eye lands on.

### 3.2 The Reference axis never renders as `danger`

§7 splits the annotation-confidence *Interpretation* column into
**Completeness** and **Reference** icons. The review specified
"success/warn/danger icons" for both.

A `danger` icon on the Reference axis for low identity would re-create
precisely the misleading signal the two-axis design exists to prevent.
Task 43b §5.4 is explicit: low identity + high coverage is *a complete
gene in an under-referenced clade* and **must not read as failure**,
and this is not hypothetical — `INT-ANIMAL-01`'s `COX3` (39.5 / 96) and
`ATP6` (39.3 / 100) are both real, intact genes. Two of three scored
genes in the fixture would show red.

**Resolved:** the Reference axis uses **success (well-referenced) and a
neutral informational marker (under-referenced) — never `danger`**. The
tooltip on the under-referenced marker carries the affirmative wording
("Complete gene in an under-referenced clade. This is not a failure."),
so the reader is *told* it is a good result rather than left to infer
it. The Completeness axis may use success/warning normally — a
truncated fragment genuinely is a caveat.

This keeps 43b §5.4 rule 3's affirmative-label requirement satisfied
through a different surface, and 43b's
`test_low_identity_high_coverage_carries_affirmative_label` should be
retargeted at the tooltip text rather than deleted.

---

## 4. Tab order — Validation moves to position 2

New strip order: **Overview, Validation, Assembly, Barcodes.**

Reasoning worth recording: trustworthiness now precedes content. A
reader who cannot trust the assembly has no reason to read its
statistics, and on a soft-failed sample Validation is already the
terminal tab (43b §3.3) — putting it immediately after Overview means
the terminal content sits where the reader arrives next.

Touch points:

- `TAB_NAMES` in `bin/report/report.py`.
- `components/content-tabs.html` currently hardcodes each pane's
  `include` against a fixed index (`tab-1` → assembly, `tab-2` →
  validation). Swap them. Consider driving the pane includes from the
  same ordered list as the buttons so the two cannot drift.
- The `warnings` mirror keys on **tab name**, not index — no change.
- Tests asserting `id="tab-2"` as Validation (§11).

**Spec update required:** spec §6a.1's tab bullet states the four tabs
answer "did it work, what did I get, should I trust it, which barcodes
can I use — *in that order*". That sentence and §6a.2's table row
numbering must be updated to the new order: *did it work → should I
trust it → what did I get → which barcodes can I use*.

---

## 5. Overview tab

1. **Layout** per §3.1: compact Input params + Run context strip at
   top, then Key findings, then warnings, then Sequencing data quality.
2. **Key findings becomes a two-column table**, not a stack of
   `alert-*` components. Column 1 the finding, column 2 the value, or a
   category/statement split — pick one and keep it uniform. The
   per-finding contextual class survives as row styling or a badge (via
   §2.1's macro), because `low_coverage` must stay visually distinct
   from `ok` without reading as a failure (43a §7).
3. **Organelle renders as a full name.** Add an `ORGANELLE_NAMES`
   constant beside the existing `ORGANELLE_TITLES` in
   `bin/report/report.py`, and use it everywhere the organelle is shown
   — currently only Input params' *Organelle* row, but the mapping
   should be the single source of truth from the outset.
   - `pt` → *Chloroplast*.
   - `mt` → the review says *Mitochondria*, which is the **plural**;
     the singular is *Mitochondrion*. One organelle genome is being
     reported, and it pairs correctly with *Chloroplast*.
     **Recommend *Mitochondrion*; confirm with the project owner
     before implementing**, as this is the owner's wording call.
4. **Status row** — no template change beyond §2.1's badge fix; the
   row is not blank, its badge is invisible.
5. **Sequencing data quality** — no code change; §2.2.

---

## 6. Assembly tab — layout

Put the per-contig table and the mean-coverage chart **side by side in
one row**:

- **Table:** minimum horizontal width to fit its content. Bootstrap's
  `.table` is `width:100%` by default, so this needs an explicit
  override (`w-auto`, plus `white-space: nowrap` on cells).
- **Chart:** **horizontal** bars (Plotly `orientation: 'h'`, contig
  names on the y-axis, coverage on the x-axis), flexing to fill the
  remaining x-space (`responsive: true`, container `width: 100%`).

Horizontal orientation is the better fit anyway: contig labels are long
and `plant_mt` runs `emit: 'all'` with `max_contigs: 20`. **Scale chart
height with contig count** rather than pinning it at 300 px, or a
20-contig sample renders unreadably compressed.

The existing bucket colouring (target / secondary / off-target) is
unchanged.

---

## 7. Assembly tab — annotation confidence behind a modal

1. **Move the whole confidence section into a modal**, triggered by an
   **"Annotation details"** button on the Assembly tab.
2. **The leading explanatory paragraph moves behind an info icon**,
   shown on hover (Bootstrap tooltip or popover — the boilerplate
   already initialises tooltips globally in `index.html`).
3. **Remove the standalone interpretation/meaning key table.**
4. **Split the *Interpretation* column into *Completeness* and
   *Reference***, each rendering an icon, each icon carrying the
   corresponding explanation as a tooltip. Severity mapping per §3.2 —
   **Reference is never `danger`.**

Mapping from the existing `bin/report/confidence.py` quadrants:

| Quadrant | Completeness | Reference |
|---|---|---|
| `complete_well_referenced` | success | success |
| `truncated_well_referenced` | warning | success |
| `complete_under_referenced` | success | info (never danger) |
| `fragmentary` | warning | info (never danger) |
| `unscored` | unscored marker | unscored marker |

**Derive the two axis labels on the result class, not in the
template** — extend `confidence.py` with the per-axis label/severity
alongside the existing quadrant label, exactly as 43b §5.4 rule 3
requires for the quadrant itself. The template renders what it is
given.

**§5.4 rule 4 is reinterpreted, not dropped.** That rule requires the
interpretation key to be visible in the tab so the reader learns the
axes. Deleting the key table satisfies the review only if the
**per-icon tooltips plus the info-icon paragraph carry that same
two-axis explanation**. Record this as a deliberate substitution; do
not delete the key and leave nothing in its place.

**Constraint carried forward:** `unscored` must remain visibly distinct
from a zero/low score (43b §5.4). An empty cell will not do — an
explicit "unscored" marker with its own tooltip is required.

Check the modal still behaves in the **saved read-only report**
(`save-report.js` clones the DOM and strips `.hide-broken` links); the
confidence table must survive that path.

---

## 8. Validation tab — Recruitment & coverage

1. **Collapse the three stacked tables into one.** Currently the gate
   table, the sibling-split table and the subsampling table render as
   separate `<table>` elements with alerts between them. One table,
   with the sibling-split and subsampling rows conditional on the same
   data that gates them today.
2. **Show the alerts above the table**, not interleaved.
3. **The sibling-organelle caveat becomes `alert-warning`**, not
   `alert-info` — this is the *"the sibling-organelle split was
   unavailable ... the estimate may over-state target depth"* message.
   **Change the severity in `build_warnings()` too**, from `info` to
   `warning`, so the Overview mirror and the Validation tab agree.
   A mirror that disagrees with the tab it mirrors is worse than
   either alone (rule 18).

The below-species-threshold note stays `info` — it flags a taxon, it
does not caveat a number.

**The merged table also gains a row, and its existing coverage row is
relabelled** — see §9, which owns both changes.

---

## 9. Validation tab — "actual coverage" (review question)

> *"If possible, the 'Estimated coverage' should be accompanied by
> 'Actual coverage' determined from RECRUIT result? (Push back on the
> user if this doesn't make sense)"*

**Pushback on the source; agreement on the substance.**

**Not from RECRUIT.** `recruit_stats.json` carries `reads_aligned`,
`reads_recruited`, `min_aligned_frac` and `min_aligned_bp` — read
*counts* and the thresholds applied. It computes no depth at any point,
so there is no "actual coverage" in it to surface.

**But the number the review is reaching for does exist, and is already
in `metadata.json`.** Flye reports a mean coverage per contig in
`assembly_info.txt`, carried through as `assembly.contigs[].coverage`.
Read against `bin_metadata.contigs_selected`, that gives observed depth
**on the contig actually selected as the target** — which is precisely
the comparison spec §2.1.6 asks the report to make, and which the
Validation tab currently states only as generic prose ("actual depth on
the assembled contig is typically well below `estimated_cov`").

All three fixtures make the point concretely:

| Sample | `estimated_cov` | Flye depth on selected contig(s) |
|---|---|---|
| `INT-ANIMAL-01` | 327.52× | 123.0× (`contig_10`) |
| `INT-PLANT-01-pt` | 254.13× | 153.0× (`contig_3`) |
| `INT-PLANT-01-mt` | 28.48× | 7.0 / 7.0 / 5.0 / 5.0× (4 contigs) |

**Add the row.** Requirements:

- **Label it honestly, and in language that needs no tool knowledge.**
  The figure is the assembler's own estimate from read assignment
  during assembly, **not** a `samtools depth` pileup — so never a bare
  *"Actual coverage"*, which would present a derived estimate as a
  measurement (rule 18).

  **Distinguish the two rows by what each was measured over**, not by
  naming the tool. The gate row above is already labelled *"Estimated
  coverage"*, so a second row reading *"Estimated coverage (Flye)"*
  gives the reader two adjacent "estimated coverage" figures whose only
  discriminator is a tool name they have no reason to recognise —
  which defeats the contrast the row exists to draw. Prefer:

  - *Estimated coverage — recruited reads* (the existing gate row)
  - *Estimated coverage — assembled mitochondrion* / *…assembled
    chloroplast* (the new row, using §5's `ORGANELLE_NAMES` mapping)

  Avoid *"mitochondrial coverage"* as the discriminator: the
  recruitment estimate is also meant to be mitochondrial coverage, so
  that phrasing implies the first number is not, which is not the
  claim. The difference is the denominator — a nominal genome size
  versus the contig actually assembled.

  **Put the tool name and the explanation in a tooltip**, roughly:
  *"Depth estimated by the assembler across the assembled organelle
  contig, rather than across the whole recruited read pool. This is
  usually lower: most recruited reads assemble into other contigs,
  which are reported as off-target on the Assembly tab."*

  Note that *"off-target reads were discarded"* would be the wrong
  gloss — they are not discarded, they are assembled into other
  contigs and shown to the reader. On `INT-ANIMAL-01` no reads are
  discarded as sibling-organelle at all (`sibling_organelle_fraction:
  0.0`), yet the estimate still falls 327× → 123×, because the
  assembly totals 215,078 bp across 10 contigs while the target contig
  is 16,799 bp of it.
- **`emit: 'all'` targets select several contigs** (`plant_mt`, up to
  20). Render a range or a per-contig list, not a single number that
  silently describes only the first.
- **`fail` and `no_assembly` have no selected contig.** Render the row
  absent or `-`, never `0`, which would read as a measured zero depth.
- **This does not close the depth-track item in `tasks/todo.md`.** That
  entry is about a *per-position* depth track needing a new alignment
  stage. This is a per-contig mean that already exists. Leave the todo
  entry in place.

The `INT-PLANT-01-mt` row is the strongest argument for the change: a
`low_coverage` sample whose estimate reads 28× while its assembled
contigs sit at 5–7×. That gap is exactly what spec §2.1.6 warns about,
and the report should show it per-sample rather than assert it in
prose.

---

## 10. Validation tab — annotation cross-checks

Same shape as §8: two tables with an alert between them. **Collapse
into one table, with the warning above it.** The genetic-code rows
(`genetic_code_annotate` / `genetic_code_cds`) and the cross-check rows
(`agreed` / `miniprot_only` / `annotator_only` /
`coordinate_conflicts`) become one table.

Preserve from 43b: canonical gene names displayed, raw
(fragment-suffixed) names exposed on hover, per task 42 §5.3.

---

## 11. Barcodes tab

The *Reason* column values render as `class="badge text-dark
bg-warning"` — i.e. via §2.1's macro with the `warning` context. This
is the one badge the review caught; the other four are §2.1.

---

## 12. Tests — extending `scripts/tests/test_report.py`

Rule 14: **100 % branch coverage** via `scripts/pytest.sh`, currently
green at 452 tests. Assert on rendered output, not on mocks.

Regression guards for §2.1 — the highest-value assertions here, because
the defect was invisible rather than loud:

- **No `text-bg-` string appears anywhere in the rendered HTML**, for
  every `sample_status`. Mechanical, cheap, and it fails the moment
  anyone reintroduces the pattern.
- **Every badge carries both a `bg-*` and a foreground class**, with
  the §2.1 pairing.
- **The sample-status badge renders its label text** for all five
  statuses.

Layout and ordering:

- **Overview DOM order** — context strip, then Key findings, then
  warnings, then read QC. 43b's existing "Key findings is inside
  `tab-0`" test must still pass unchanged (§3.1).
- **Tab order** — `Validation` is `tab-1` and `Assembly` is `tab-2`;
  Overview is still `tab-0` and active on load. Update 43b's
  `test_low_coverage_renders_warn_floor_warning_in_both_places`, which
  currently reads Validation out of `id="tab-2"`.
- **Key findings renders as a table**, not `alert-*` list items, and
  `low_coverage` remains visually distinct from `ok` without failure
  wording (43a §7 — assert the existing guarantee survives the
  refactor).

Content:

- **`ORGANELLE_NAMES` mapping** — `mt` and `pt` each render their full
  name; an unrecognised organelle degrades gracefully rather than
  raising.
- **Reference axis is never `danger`** — assert directly against
  `COX3` (39.5 / 96) and `ATP6` (39.3 / 100) from `INT-ANIMAL-01`.
  This is the §3.2 regression and the one that matters most.
- **The affirmative "under-referenced" wording is present** in the
  tooltip after the key table is removed — retarget 43b's existing
  assertion rather than dropping it.
- **`unscored` renders its own marker**, distinct from any zero or
  low-score rendering.
- **The quadrant → (completeness, reference) axis mapping is
  unit-tested directly** on `confidence.py`, independent of rendering,
  for all five quadrants.
- **Annotation details modal** is present and contains the confidence
  table.
- **Sibling caveat is `warning` severity in both** the Validation tab
  and the `build_warnings()` mirror — assert the two agree.
- **Flye depth row** (§9): present and labelled for a single-contig
  target; renders every contig for an `emit: 'all'` target; absent or
  `-` (never `0`) for `fail` / `no_assembly`.

Unchanged and must stay green: self-containment (no external fetch)
across every status, with 43b's two allowlisted strings.

---

## 13. Sequencing relative to task 45

**This task should be executed before task 45_run_report.md.** Task 45
reuses the shared report components, and §2.1's badge macro plus §5's
organelle-name mapping are exactly the kind of thing it would otherwise
fork in its broken form — shipping the invisible-badge defect into the
run-level report as well.

**Renumbering was considered and rejected.** Task 45 is referenced by
name in spec 06a-reports.md and in the completed briefs for tasks 43a
and 43b; renaming it to keep numeric order would stale those
references for no functional gain. The dependency is recorded here
instead.

---

## 14. Spec updates

- **spec/06a-reports.md §6a.1** — the four-tab bullet's question order
  becomes *did it work → should I trust it → what did I get → which
  barcodes can I use* (§4). §6a.2's table rows renumber to match.
- **spec/06a-reports.md §6a.2** — the Overview row's block ordering
  gains the context strip ahead of Key findings (§3.1). Key findings
  remains described as above the fold; do not weaken that wording.
- **No change to §6a.5** — the renderer's structure is unchanged.
- **tasks/todo.md** — unchanged. In particular, **do not delete the
  per-position depth-track entry** on the strength of §9's per-contig
  mean (§9, final bullet).

Per the task-authoring convention: reference tasks in spec as plain
text (`task 46_report_enhancement.md`), never as markdown links, and
never with line anchors.

---

## 15. Acceptance criteria

1. No `text-bg-*` class appears in any template or rendered report;
   all five badges render legibly via a single shared macro, asserted
   by test (§2.1).
2. The Overview tab renders context strip → Key findings → warnings →
   read QC, with Key findings above the fold and 43b's DOM-order test
   still passing (§3.1).
3. Key findings renders as a two-column table, with `low_coverage`
   still distinguishable from `ok` and still not framed as a failure.
4. The tab strip reads Overview, Validation, Assembly, Barcodes;
   Overview active on load.
5. Organelle renders as a full name from one shared mapping, with the
   `mt` wording confirmed by the project owner (§5 item 3).
6. The Assembly tab's contig table and horizontal coverage chart sit
   side by side, the chart flexing to available width and scaling its
   height with contig count.
7. Annotation confidence lives in an "Annotation details" modal, with
   Completeness and Reference as separate icon columns; **the Reference
   axis never renders `danger`**, and the affirmative under-referenced
   wording is present in a tooltip (§3.2).
8. Recruitment & coverage and Annotation cross-checks each render as
   one table with alerts above; the sibling caveat is `warning`
   severity in both the tab and the Overview mirror.
9. Flye mean depth on the selected target contig(s) renders alongside
   `estimated_cov`, honestly labelled, handling `emit: 'all'` and the
   no-contig statuses (§9).
10. `flake8` clean at 79 columns; `scripts/pytest.sh` green at 100 %
    branch coverage; report remains a single self-contained file.
11. `-stub-run` green; `tests/integration/assertions.sh` still passing.
12. A fresh `-profile integration` run is performed and the report
    re-reviewed against non-stale output, confirming `read_qc` and
    `coverage.recruitment` populate (§2.2).

---

## 16. Out of scope

- **Upgrading the vendored Bootstrap** — §2.1; the macro is the fix.
- **A per-position coverage depth track** — still a new alignment
  stage and a new task; §9 does not close it.
- **The real organelle map SVG** — task 44_organelle_map.md.
- **`run-report.html`** — task 45_run_report.md, which this task
  should precede (§13).
- **Re-sourcing fixtures for `fail` / `no_assembly` / `no_barcode`** —
  the known gap carried since task 42 is unchanged here.

---

## 17. Outcomes

- **§5 item 3 (mt wording) confirmed with the project owner before
  implementation**, as the brief required: **"Mitochondrion"**
  (singular), not the review's "Mitochondria". `ORGANELLE_NAMES` in
  `bin/report/report.py` is the single source of truth (`mt` →
  Mitochondrion, `pt` → Chloroplast), consumed by the Input params
  row and the new Validation-tab Flye-depth row label.
- **§2.1 badge defect.** `scripts/report/templates/macros/badge.html`
  (`render_badge(text, context, title='')`) owns the
  `bg-*`/foreground pairing from the task's luminance table and is
  now the only place a badge's classes are constructed. All five call
  sites converted (`inputs.html` sample status, `validation.html`
  gate decision + below-species-threshold + annotator-only
  cross-check, `barcodes.html` drop-out reason). `grep -c
  'text-bg-'` across `scripts/report/templates/` returns 0 (excluding
  the macro's own explanatory comment).
- **§2.2 stale fixture diagnosis confirmed independently.** Running
  `tests/integration/assertions.sh` against the existing (gitignored,
  pre-43a) `tests/integration/output/` showed all three sample
  `report.html` files are **zero bytes** — they predate task 43a's
  real renderer entirely, consistent with the brief's account. This
  is pre-existing state, not something introduced or fixed by this
  task; a fresh `-profile integration` run (§12 acceptance criterion,
  deferred — see below) is what regenerates them.
- **§4 tab reorder.** `TAB_NAMES` and a new `TAB_TEMPLATES` mapping in
  `report.py` drive both the button strip and the pane markup in
  `content-tabs.html` from one ordered loop, so the two structures
  cannot drift the way the old hardcoded `tab-1`/`tab-2` includes
  could.
- **§3.1/§5 Overview relayout.** `overview.html` now orders context
  strip → Key findings → warnings mirror → read QC. `inputs.html`,
  `walltime.html` and `provenance.html` all use the existing
  `tight lined font-small` styling so the strip reads as visually
  subordinate to Key findings.
- **§5.2 Key findings as a table.** Each finding in `key_findings()`
  gained a `label` field (Outcome/Assembly/Coverage/Annotation/
  Barcodes/Top hit); `key-findings.html` renders a two-column table
  with `table-{{ finding.class }}` row styling instead of stacked
  `alert-*` list items.
- **§6/§7 Assembly tab.** Per-contig table (`w-auto`, `white-space:
  nowrap`) and the coverage chart now share one row; the chart is
  horizontal (`orientation: 'h'`) with height `Math.max(300, n *
  30)`. Annotation confidence moved into an "Annotation details"
  modal; the explanatory paragraph moved behind an info-icon tooltip;
  the standalone interpretation key table was deleted. `confidence.py`
  gained `axis_severity()`/`axis_tooltip()` and an `AXIS_SEVERITY`
  table implementing §7's quadrant mapping — the Reference axis is
  `success`/`info`/`unscored` only, never `danger`. A new
  `macros/confidence-icon.html` renders each axis as an emoji icon
  with its tooltip carrying the two-axis explanation the deleted key
  table used to.
- **§8/§10 Validation tab collapses.** The gate/sibling-split/
  subsampling tables merged into one table with both alerts moved
  above it; the annotation cross-check tables merged the same way.
  The sibling-organelle caveat is now `warning` severity in both
  `build_warnings()` and the inline Validation-tab alert (previously
  `info`), so the Overview mirror and the tab agree.
- **§9 Flye depth row.** `validation_view()` gained `flye_depth`
  (`_flye_depth()`), reading `bin_metadata.contigs_selected` against
  `assembly.contigs[].coverage`. Row label is *"Estimated coverage —
  assembled {mitochondrion|chloroplast}"* with the tool-name/
  explanation in a tooltip; renders every contig for `emit: 'all'`
  targets (comma-joined); absent entirely (not `0`) when no contig
  was selected (`fail`/`no_assembly`).
- **§11 Barcodes tab** — the one badge the review caught now also
  goes through the shared macro.
- **§12 tests.** `scripts/tests/test_report.py` grew from 452 to
  **472 cases** (11 new test classes covering the badge regression,
  Overview layout order, Key-findings-as-table, tab order, organelle
  names, the confidence axis mapping including the real
  COX3/ATP6 `INT-ANIMAL-01` values, the annotation-details modal, the
  sibling-caveat severity agreement, and the Flye depth row's three
  cases). Two 43b tests were retargeted rather than deleted, per the
  brief: `test_low_identity_high_coverage_carries_affirmative_label`
  now asserts the affirmative wording in the Reference-axis tooltip
  instead of a visible key-table row, and
  `test_low_coverage_renders_warn_floor_warning_in_both_places` reads
  Validation out of `tab-1` instead of `tab-2`. `scripts/pytest.sh
  scripts/tests/` — **472 passed, 100% branch coverage** across every
  `bin/*.py` and `bin/report/**` module. `flake8` clean at 79 columns.
- **§14 spec updates.** `spec/06a-reports.md` §6a.1's tab-order
  sentence and §6a.2's table (row order and the Overview row's block
  description) updated to match. §6a.2's Assembly row description
  updated for the modal/two-axis-icon change. No change to §6a.5, as
  specified.
- **§13 sequencing.** This task lands before task 45_run_report.md, as
  required — the badge macro and `ORGANELLE_NAMES` are now the shared
  components task 45 will reuse rather than fork.
- **Deviations from the brief / deferred work:**
  - **A fresh `-profile integration` run (acceptance criterion 12) was
    deferred**, by explicit user decision during execution (the run
    is real-tool/network/Docker-heavy and was judged out of proportion
    to this task's template/render-logic scope). `-stub-run -resume`
    was run instead and completed green end-to-end (exit 0, including
    `COLLATE` and `RUN_REPORT`), and the real `INT-ANIMAL-01`
    `metadata.json` fixture was rendered directly through the updated
    templates to confirm no exceptions and that the new content
    (organelle name, Flye depth row, affirmative tooltip, confidence
    modal, zero `text-bg-` occurrences) all appear correctly. The
    live re-run that would regenerate the stale `report.html` fixtures
    and confirm `read_qc`/`coverage.recruitment` populate on fresh
    output remains outstanding — carried into `tasks/todo.md`.
  - No other deviations from the brief.
