# Task 53 — Bootstrap-styled hover tooltips for the assembly graph

**Phase:** P4a (spec 06-phases.md).
**Goal:** Replace the assembly graph's plain native browser tooltip
(an SVG `<title>` child, task 47_bandage_interactive_svg.md) with the
same Bootstrap tooltip styling the rest of `report.html` already uses
(the info badges, the confidence icons), without giving up the
no-JavaScript fallback task 47 deliberately built in.

**Depends on:** task 47_bandage_interactive_svg.md (ships the
`data-node`/`data-contigs`/`<title>` markup this task upgrades).

**Related tasks:** task 52_graph_bucket_colouring.md. Independent —
either can land first. Task 52 owns the `<title>` *text* (it appends a
bucket suffix such as `(target)` / `(mixed: target, off-target)` in
C12); this task re-presents whatever text C12 emits, so the two
compose with no extra work. The only shared file is
`scripts/report/templates/components/assembly.html` (task 52 edits the
graph info badge and adds a colour key; this task adds one script
call), so whichever lands second should re-read the other's diff to
that file first.

---

## 1. Overview

`report.html` already has a consistent hover-tooltip system: the info
badges (`(i)` icons next to section headings) and the
annotation-confidence icons both use Bootstrap's tooltip component — a
small dark popover, themed with the rest of the page and positioned
around the trigger element (see
`scripts/report/templates/macros/info-badge.html`,
`macros/confidence-icon.html`, and the global initialiser at the
bottom of `scripts/report/templates/index.html`).

The assembly graph's node tooltips (task 47) don't use this. They rely
on the SVG `<title>` child element — the plain browser tooltip every
SVG viewer supports natively, with no JavaScript involved. That's
slow to appear, unstyled, and looks different from every other
tooltip on the page. It was a deliberate choice for the *published
file*, though: `graph.svg` is also copied verbatim into
`diagnostics/graph.svg` as a standalone diagnostic, and task 47's
acceptance criterion 3 required the hover label to work there with no
JavaScript. That's a requirement on the SVG file being self-contained,
not a rule against the *report page* layering something nicer on top
when it already has Bootstrap loaded.

This task is that layering: inside `report.html`, graph node tooltips
look and behave like every other tooltip on the page. Opened on its
own, the same `graph.svg` keeps working exactly as task 47 left it.

---

## 2. What the current code already gives this task

Verified against the current tree:

- **Bootstrap is v5.0.2** (`scripts/report/static/js/bootstrap.bundle.min.js`).
  This matters, because several tooltip features people reach for
  arrived later:
  - Content comes only from the trigger's `title` **attribute** (or a
    `title` config option). `data-bs-title` is **not** supported in
    5.0.2 (it arrived in 5.2) — don't use it.
  - On construction, `_fixTitle()` moves `title` into
    `data-bs-original-title`, sets `title=""` (so the browser shows
    nothing natively), and sets `aria-label` from it when the element
    has no text content. Content is read back from
    `data-bs-original-title` at show time.
  - The tooltip auto-hides on `hide.bs.modal` **only for the modal its
    element was inside when the tooltip was constructed**
    (`this._element.closest('.modal')`, bound once). See §3.4 for why
    that matters here.
  - `bootstrap.Tooltip.getInstance(el)` exists, and so does `.hide()`.
  - Default `container` is `document.body`; tooltip z-index (1080 in
    the bundled `bootstrap.min.css`) sits above the modal (1060), so a tooltip over the fullscreen
    graph modal is visible.
- `index.html`'s global initialiser (after `</body>`) runs once at
  load: `querySelectorAll('[data-bs-toggle="tooltip"]')` → `new
  bootstrap.Tooltip(el)`. Every `components/*.html` template, including
  `assembly.html`, is parsed before it, so anything an inline script in
  `assembly.html` does to the graph's DOM is in place when the scan
  runs.
- `bin/report/report.py`'s `get_static_file_contents()` inlines every
  file under `scripts/report/static/js/` into `<head>`.
  `modal-move-restore.js` (task 47) is the precedent for a small
  shared helper defined there and called from an inline `<script>` in
  a component template.
- `bin/annotate_graph_svg.py` (task 47) puts the `<title>` as the
  **first child** of each `<g data-node="..." data-contigs="...">`.
- The graph thumbnail `<div id="assembly-graph">` carries
  `data-bs-toggle="modal"` and opens `#graphModal`. Its children (the
  `<svg>`) are *moved* into `#assembly-graph-modal-body` on
  `show.bs.modal` and back on `hidden.bs.modal`
  (`bindSvgModalMoveRestore`). Node `<g>` elements with their own
  `data-bs-toggle="tooltip"` are nested inside that div. Bootstrap's
  modal click handler is delegated and walks up to the nearest
  `[data-bs-toggle="modal"]`, so clicking a node still opens the
  modal — the inner attribute doesn't shadow it.
- `scripts/report/static/js/save-report.js` saves the page by cloning
  the **live DOM** (`document.documentElement.cloneNode(true)`). So a
  saved report contains the post-conversion, post-Bootstrap markup
  (`title=""`, `data-bs-original-title`, `data-bs-toggle="tooltip"`,
  no `<title>` child), and must still work when reopened (§3.5).
- The organelle map (`bin/organelle_map.py`, task 44) also uses
  native `<title>` children, inside `#organelle-map`, with the same
  modal move/restore pattern.

---

## 3. Design

### 3.1 Keep C12 itself unchanged

`bin/annotate_graph_svg.py` must **not** learn about Bootstrap or
`data-bs-toggle`. Its output is also published as a standalone
diagnostic (`diagnostics/graph.svg`), where that markup would mean
nothing. C12's contract stays as task 47 (and, if landed, task 52) left
it: `data-node`, `data-contigs`, a plain `<title>` child. This task
lives entirely in the report's static JS + template layer.

### 3.2 A shared helper: promote `<title>` children to Bootstrap tooltips

New file `scripts/report/static/js/svg-title-tooltips.js`, one
function, following `modal-move-restore.js`'s style (header comment
explaining why it exists; plain ES5-ish JS like its neighbour):

```
promoteSvgTitleTooltips(containerId, itemSelector):
    container = getElementById(containerId); return if missing
    for each el matching itemSelector under container:
        titleEl = el's direct <title> child (not a descendant's)
        if titleEl:
            el.setAttribute('title', titleEl.textContent)
            titleEl.remove()
            el.setAttribute('data-bs-toggle', 'tooltip')
```

Called from an inline `<script>` in `assembly.html`, placed after the
`#assembly-graph` div (so the SVG is parsed) and before the page's
global init (always true — see §2):

```
promoteSvgTitleTooltips('assembly-graph', 'g[data-node]')
```

Why this shape:

- **No double tooltip.** Removing the `<title>` child removes the
  native tooltip; Bootstrap's `_fixTitle()` then blanks the `title`
  attribute too. There's no window where both show. (Setting
  `data-bs-toggle` alongside a surviving `<title>` child would leave
  the native one active.)
- **No new init path.** Setting `data-bs-toggle="tooltip"` before the
  global scan runs means the existing initialiser in `index.html`
  constructs the tooltips. Nothing changes in `index.html`.
- **Direct child only** — so a future SVG with nested titled elements
  can't have an inner `<title>` hoisted onto an outer group.
- **Scoped** by container id + selector, so the organelle map can
  adopt it later with one line (§3.6) and nothing else on the page is
  touched.
- **Idempotent** on a saved report: the `<title>` children are already
  gone, so it's a no-op (§3.5).

Use `data-bs-placement="top"` (set in the same loop) to match the info
badges, unless it clips at the top edge of the fullscreen modal — then
let Popper flip (the default fallback behaviour already does this).

### 3.3 Tooltip content

Use exactly the `<title>` text C12 built (segment name, or
`segment — contig,contig`, plus task 52's bucket suffix if landed).
Don't rebuild it client-side from `data-node`/`data-contigs` — one
source of truth for the wording, and the standalone file and the
report say the same thing.

Tooltip text goes through the `title` attribute as plain text.
`html: false` (the default) must stay, so a contig name can never be
interpreted as markup.

### 3.4 Modal move/restore: hide tooltips explicitly around the move

This is the one real wrinkle. Tooltip instances are tied to the
element, not its position in the DOM, so they **keep working after
the SVG moves into the modal and back** — no dispose/re-create is
needed for the tooltips to function. But two stale-tooltip cases
exist, both from §2's 5.0.2 behaviour:

1. **Thumbnail → modal.** The user hovers a node (tooltip showing) and
   clicks it to open the modal. The node moves into the modal body,
   but its tooltip stays shown, anchored to where the node *was*,
   until the pointer moves.
2. **Modal → thumbnail.** The user hovers a node in the fullscreen
   modal and closes it (Esc / close button). Bootstrap's built-in
   "hide tooltips when their modal hides" hook was bound at
   construction, when the node was in the thumbnail with no `.modal`
   ancestor — so it never fires. The tooltip can be left floating
   over the page after the modal is gone.

Fix: in the same helper file, a second small function — or an
optional argument to `promoteSvgTitleTooltips` — that, given the modal
id, on both `show.bs.modal` and `hide.bs.modal`, calls
`bootstrap.Tooltip.getInstance(el)` → `.hide()` for every promoted
element in the container. Hook `hide.bs.modal` (before the fade), not
`hidden.bs.modal`, so nothing lingers during the transition.
`modal-move-restore.js` itself stays unchanged — it moves DOM nodes
and shouldn't know about tooltips.

**Verify both cases by hand** (§5.1) rather than assuming. If
Popper positioning turns out wrong inside the fullscreen modal after
the move (it shouldn't — Popper recomputes on show), call `.update()`
on the instance in `shown.bs.modal` as a fallback. If that still
fails, dispose and re-create the instances around the move.

### 3.5 Saved reports

Because `save-report.js` clones the live DOM, a saved report contains
Bootstrap's post-construction attributes and no `<title>` children.
On reopening:

- `promoteSvgTitleTooltips` finds no `<title>` children → no-op.
- The global scan still finds `data-bs-toggle="tooltip"`;
  `_fixTitle()` sees `title=""` with an existing
  `data-bs-original-title` and leaves it alone; the tooltip shows the
  original text.

This is the same path the info badges already go through on save, so
it should just work — but it's part of the manual check (§5.1) because
this is the first time it's applied to SVG elements. Also note: if a
tooltip is *open* when the user saves, the clone carries its
`aria-describedby` and the tooltip `<div>` in `<body>`. That's
pre-existing behaviour shared with the info badges, and out of scope
here — don't fix it in this task, but mention it in Outcomes if seen.

### 3.6 Scope

Only the assembly graph (`#assembly-graph`). The organelle map
(`#organelle-map`, task 44_organelle_map.md) has the same
native-`<title>`-only tooltips and would gain consistency from the
same one-line call plus the §3.4 modal hook. It's deliberately left
out to keep this task to what was asked. Add it to `tasks/todo.md` as
a one-line follow-up ("apply `promoteSvgTitleTooltips` to
`#organelle-map` — check its `<title>` placement first").

### 3.7 Accessibility note

After conversion each node gets an `aria-label` (from `_fixTitle`),
which is at least as good as the native `<title>` for screen readers.
The nodes are not keyboard-focusable (`<g>` without `tabindex`), so
Bootstrap's `focus` trigger does nothing for them. That's unchanged
from task 47 and stays out of scope (§6).

---

## 4. Work items

1. `scripts/report/static/js/svg-title-tooltips.js` (new): the
   §3.2 conversion and the §3.4 modal hide hook.
2. `scripts/report/templates/components/assembly.html`: one inline
   `<script>` call for `#assembly-graph` / `#graphModal`, next to the
   existing `bindSvgModalMoveRestore('graphModal', ...)` call. The
   conversion must run whether or not the modal exists, so if it's
   placed inside the `{% if graph_svg %}` modal block, confirm that
   block is always rendered when the graph div is.
3. Confirm `index.html`'s global init picks the nodes up with no
   change there.
4. `tasks/todo.md`: the organelle-map follow-up line (§3.6).
5. Manual verification (§5.1), recorded in this task's Outcomes.

---

## 5. Tests

No `bin/*.py` changes, so no new `scripts/pytest.sh` coverage is
required for Python.

### 5.1 Manual verification (required — record results in Outcomes)

Against a report rendered from a real integration sample with a
multi-node graph (`INT-PLANT-01-pt` is the best candidate — the IR
gives several nodes). Check in at least Chrome and Firefox:

1. Hover a thumbnail node → Bootstrap-styled tooltip with the C12
   text; no native tooltip appears, even after holding still ~2 s.
2. Hover a node, click it → modal opens, no tooltip left behind on the
   page.
3. In the modal, hover nodes → tooltips correctly positioned over the
   enlarged graph and visible above the modal.
4. Hover a node in the modal, press Esc → no tooltip left floating.
5. After closing, hover thumbnail nodes again → still working. Repeat
   open/close twice more.
6. Save the report (`save-report.js`), open the saved file, repeat 1–5.
7. Open `diagnostics/graph.svg` directly in a browser → native
   `<title>` tooltip still works.

### 5.2 Unit — `scripts/tests/test_report.py` (small additions)

- The rendered page inlines `svg-title-tooltips.js` and calls
  `promoteSvgTitleTooltips` when a graph SVG is present, and doesn't
  call it when `graph_svg` is absent (template branch coverage).
- The graph markup embedded in `report.html` still contains the
  `<title>` children at render time — the conversion is client-side,
  so the server-rendered HTML must be unchanged from task 47.
- Existing self-contained assertions
  (`test_self_contained_for_every_status` and friends) keep passing —
  no external asset reference introduced.

### 5.3 Integration — `tests/integration/assertions.sh`

No changes. The `graph.svg` assertion reads
`diagnostics/graph.svg`, which this task leaves byte-for-byte
unchanged (§3.1), so it must keep passing as is. A fresh integration
run isn't required for this task alone. Re-rendering the report from
existing outputs is enough for §5.1. Give the rendered report path to
the user (CLAUDE.md).

---

## 6. Explicit non-goals

- Any change to `bin/annotate_graph_svg.py` or the published
  `diagnostics/graph.svg` (§3.1).
- The organelle map's tooltips (§3.6) — follow-up only.
- Upgrading Bootstrap to get `data-bs-title` or other ≥5.2 tooltip
  features — work within 5.0.2.
- Interaction beyond hover: click-to-pin, keyboard focus on nodes,
  node highlighting — same boundary task 47 §9 drew.
- Fixing the pre-existing "open tooltip captured by save" behaviour
  (§3.5).

---

## 7. Acceptance criteria

1. Hovering a graph node inside `report.html` shows a Bootstrap-styled
   tooltip, the same look as the info badges and confidence icons,
   carrying exactly C12's `<title>` text.
2. No native browser tooltip appears on graph nodes in the report at
   any point (no double tooltip).
3. No tooltip is left behind when the fullscreen graph modal opens or
   closes. Tooltips keep working in both the thumbnail and the modal
   across repeated open/close cycles.
4. A report saved via the Save button and reopened behaves the same as
   1–3.
5. `diagnostics/graph.svg` is unchanged by this task and still shows
   its native `<title>` tooltip when opened on its own.
6. `report.html` remains fully self-contained. `test_report.py`
   passes, including the §5.2 additions; `assertions.sh` passes
   unchanged.
7. Manual verification results (§5.1) are recorded in this task's
   Outcomes section.
