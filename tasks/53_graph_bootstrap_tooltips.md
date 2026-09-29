# Task 53 — Bootstrap-styled hover tooltips for the assembly graph

**Phase:** P4a (spec 06-phases.md).
**Goal:** Replace the assembly graph's plain native browser tooltip
(an SVG `<title>` child, task 47_bandage_interactive_svg.md) with the
same Bootstrap tooltip styling the rest of `report.html` already uses
(the info badges, the confidence icons), without giving up the
no-JavaScript fallback task 47 deliberately built in.

**Depends on:** task 47_bandage_interactive_svg.md (ships the
`data-node`/`data-contigs`/`<title>` markup this task upgrades).
Independent of task 52_graph_bucket_colouring.md — either can land
first; if task 52 lands first, its bucket name can be folded into the
tooltip text (§3.3), otherwise this task ships without it and task 52
extends the text later.

---

## 1. Overview

`report.html` already has a working, consistent hover-tooltip system:
the info badges (`(i)` icons next to section headings) and the
annotation-confidence icons both use Bootstrap's tooltip component —
a small styled popover, themed to match the rest of the page, positioned
intelligently around the trigger element (see
`scripts/report/templates/macros/info-badge.html` and
`macros/confidence-icon.html`, and the global initialiser at the
bottom of `scripts/report/templates/index.html`).

The assembly graph's node tooltips (task 47) don't use this. They rely
on the SVG `<title>` child element — the plain, unstyled browser
tooltip every SVG viewer supports natively, with no JavaScript
involved. That was a deliberate choice for the *published artefact*:
`graph.svg` is copied verbatim into `diagnostics/graph.svg` as a
standalone file, and task 47's acceptance criterion 3 required the
hover label to work with **no JavaScript involved in producing it** —
a requirement about the SVG file being self-contained and inspectable
on its own, outside the report page, not a rule against the *report
page* layering something nicer on top when it has Bootstrap loaded
anyway.

This task is that layering: when `graph.svg` is viewed inside
`report.html`, its tooltips should look and behave like every other
tooltip on the page. When the same file is opened on its own (from
`diagnostics/graph.svg`, or pasted into another document), it should
keep working exactly as task 47 left it — native, plain, no JS
required.

---

## 2. What the current code already gives this task

- `scripts/report/templates/index.html` already scans the whole
  rendered page once, at load, for `[data-bs-toggle="tooltip"]` and
  initialises a `bootstrap.Tooltip` on each — the same mechanism the
  info badges and confidence icons already use. No new global wiring
  needed, just getting the graph's nodes into that scan correctly.
- `bin/annotate_graph_svg.py` (task 47) already attaches a `<title>`
  child to every annotated `<g data-node="..." data-contigs="...">` —
  the exact text this task needs to source the Bootstrap tooltip's
  content from.
- Bootstrap's tooltip reads its content from the trigger element's
  `title` HTML attribute (or `data-bs-title`) at construction time —
  **not** an SVG `<title>` child. The two mechanisms don't know about
  each other, which is the wrinkle this task has to resolve (§3.2).

---

## 3. Design

### 3.1 Keep C12 itself unchanged

`bin/annotate_graph_svg.py` should **not** be taught about Bootstrap
or `data-bs-toggle` — that would leak a report-presentation concern
into a component whose output is also published as a standalone
diagnostic file (`diagnostics/graph.svg`), and would mean the *emitted
SVG* carries markup that means nothing outside the report page. Keep
C12's contract exactly as task 47 left it: `data-node`, `data-contigs`,
a plain `<title>` child, nothing else. This task's change lives
entirely in the report template/JS layer, which is free to add
presentation on top of a plain artefact it doesn't own the generation
of.

### 3.2 Converting `<title>` to a Bootstrap tooltip without a double tooltip

Simply adding `data-bs-toggle="tooltip"` next to the existing `<title>`
child would risk two overlapping tooltips (the browser's native one
from `<title>`, and Bootstrap's own) — undefined-feeling and
inconsistent across browsers. Resolve it client-side, once, before
Bootstrap's own init scan runs: for each annotated graph node, read the
`<title>` child's text, move it onto a `title` attribute, remove the
`<title>` child, then let the existing global scan in `index.html`
pick the element up as usual (it already looks for
`data-bs-toggle="tooltip"`, so this step also needs to set that
attribute). Pseudocode, to live in
`scripts/report/templates/components/assembly.html` alongside the
existing modal move/restore script, and to run *before* the page's
existing tooltip-init scan (script tags execute in document order, and
this one is already earlier in the page than that scan):

```
for each <g data-node> under #assembly-graph:
    title_el = its direct <title> child
    if title_el exists:
        set attribute title = title_el.textContent
        remove title_el
        set attribute data-bs-toggle = "tooltip"
```

This only ever runs once, on page load, against the graph as rendered
server-side — it never needs to re-run when the modal move/restore
script (task 47 §4.6) relocates the same DOM nodes into the fullscreen
modal and back, since Bootstrap's tooltip instance is associated with
the element itself, not its position in the tree. **Verify this
assumption during implementation** (open the fullscreen graph modal,
hover a node, close it, hover the thumbnail version again) rather than
taking it on faith — if a moved node's tooltip stops working, the
instance needs disposing and re-creating around the move, not just the
DOM node relocating.

### 3.3 Tooltip content

Reuse exactly the text task 47's `<title>` already builds (segment
name, or `segment — contig,contig` when known) — do not re-derive it
differently client-side. If task 52_graph_bucket_colouring.md has
landed by the time this task starts, its bucket-colour meaning can be
folded into the same tooltip text (e.g. appending "(target)" /
"(mixed: target, off-target)") so the colour and the tooltip agree;
if not, ship without it — this task does not depend on task 52 and
should not block on it.

### 3.4 Scope

Limited to the assembly graph (`#assembly-graph` in
`scripts/report/templates/components/assembly.html`), matching what
was actually asked for. The organelle map (task 44_organelle_map.md)
has the same native-`<title>`-only tooltip on its own features and
would visually benefit from the same treatment for consistency, but
that is a separate, smaller follow-up, not required here — note it as
a one-line "also worth doing" rather than pulling it into this task's
scope.

---

## 4. Work items

1. Small script in `scripts/report/templates/components/assembly.html`
   implementing §3.2's conversion, scoped to `#assembly-graph`'s
   annotated nodes only.
2. Manual verification of the modal move/restore interaction (§3.2)
   with a real rendered report, both directions (thumbnail → modal →
   thumbnail), documented in this task's Outcomes.
3. Confirm the existing global tooltip-init scan in
   `scripts/report/templates/index.html` picks the converted nodes up
   with no changes needed there.

---

## 5. Tests

This task is report-template/JS only — no `bin/*.py` changes, so no
new `scripts/pytest.sh` coverage is required. Verification is manual
against a rendered report (§4.2) plus the existing automated checks
that must keep passing unchanged:

- `scripts/tests/test_report.py` — no new cases required, but the
  existing self-contained-file assertions
  (`test_self_contained_for_every_status` and friends) must keep
  passing: this task must not introduce any external asset reference.
- `tests/integration/assertions.sh` — the existing `graph.svg`
  structural assertion (task 47 §5.3) reads `diagnostics/graph.svg`,
  the standalone copy, which this task must leave byte-for-byte
  unchanged (§3.1) — it should keep passing with no changes.

---

## 6. Explicit non-goals

- Any change to `bin/annotate_graph_svg.py` or the published
  `diagnostics/graph.svg` file (§3.1) — this task is report-rendering
  only.
- The organelle map's tooltips (§3.4) — noted as a follow-up, not done
  here.
- Any interaction beyond hover (click-to-pin, keyboard focus styling
  beyond whatever Bootstrap's tooltip already provides) — out of
  scope, same boundary task 47 §9 already drew.

---

## 7. Acceptance criteria

1. Hovering a graph node inside `report.html` shows a Bootstrap-styled
   tooltip (consistent look with the info badges / confidence icons),
   not the plain native browser tooltip.
2. `diagnostics/graph.svg` (the standalone published file) is
   unchanged by this task — still a plain `<title>`-only tooltip, no
   Bootstrap markup, still opens and shows its native tooltip when
   viewed outside the report.
3. No double tooltip (native + Bootstrap) is visible at any point.
4. The tooltip keeps working after opening and closing the graph's
   fullscreen modal, in both the thumbnail and modal contexts.
5. `report.html` remains fully self-contained — no new external asset
   references (existing `test_self_contained_for_every_status`
   assertions keep passing).
6. `tests/integration/assertions.sh` stays fully green with no changes
   needed to the `graph.svg` assertion.
