# Task 43b — Per-sample `report.html`: restructure the tab strip and fill the panes

**Phase:** P4a (spec 06-phases.md).
**Goal:** Replace 43a's eight-stage tab strip with spec 06a-reports.md
§6a.2's **four reader-perspective tabs** — Overview, Assembly,
Validation, Barcodes — and fill them.

**Depends on:** task 43a_report_scaffold.md. It fixed the renderer's CLI
surface, the `warnings` mirror mechanism, the Key findings generator and
every plumbing input including `metadata.json`'s `read_qc` block.

**Related tasks:** 44_organelle_map.md supplies the real SVG for the
Assembly tab; this task builds against its zero-byte stub and inlines
whatever it is given, so task 44's tool decision is not a blocker.
45_run_report.md reuses the shared components; factor rather than inline.

---

## 1. Overview

43a answered "did this sample work". This task answers "what exactly do
I have", "should I trust it" and "which barcodes can I use".

Almost all of it is presentation of data already in `metadata.json`.
Very little new logic is needed. What makes this task worth a brief of
its own is §5: eight specific ways a reader can be misled by a correct
number rendered carelessly, each traceable to a CONSTITUTION principle
and each individually testable. That section is the task; the tab layout
is the easy part.

---

## 2. The strip changes shape — this is not just filling panes

**43a registered eight tabs, one per pipeline stage, in pipeline order.
That layout is superseded** (spec §6a.1, revised 2026-09-03 with the
project owner). The pipeline's stage boundaries are an implementation
detail; a biosecurity officer asks four questions, and the tabs answer
those, in that order. Several stages therefore share a tab, and the
coverage gate is deliberately filed under the question it answers
(*should I trust this*) rather than the stage that produced it.

43a's §12 acceptance criterion 5 ("the eight-tab strip is registered in
pipeline order") and criterion 7 ("the renderer's complete CLI surface
is fixed — the contract does not change shape twice") are both
superseded by this task. **Record that plainly** in 43a's completed
outcomes rather than quietly diverging from it: the contract *is*
changing shape a second time, and the reason is a design decision made
after seeing the rendered output, not a defect in 43a's execution.

### 2.1 Scaffold corrections carried into this task

Four changes to what 43a shipped, all agreed 2026-09-03:

1. **Organelle-specific H1 title.** `organelle: "mt"` → *Mitochondrial
   genome assembly*; `organelle: "pt"` → *Chloroplast genome assembly*.
   Mapping is a constant in `bin/report/report.py`, not template logic.
   Replaces the current fixed "Organelle assembly + barcode recovery
   report".
2. **Everything below the heading block moves into a tab.** 43a put the
   Key findings / Input parameters two-column row above the strip,
   always visible. It now sits at the top of the Overview tab. The
   heading block (H1, sample ID, logo, subtitle) is the only
   always-visible content.
3. **Restore the wall-time component.** 43a deleted
   `templates/components/walltime.html` because its
   `*_start_timestamp.txt` glob had no source in this repo (43a §14).
   It comes back, with a real source: `workflow.start` passed through
   `COLLATE` as a CLI argument, end time taken at render, wall time
   derived. Facility and analyst come from two new params (§2.2).
4. **`barcodes.fasta` becomes a renderer argument.** The Barcodes tab
   offers the recovered barcodes as a FASTA download; `collate.py`
   already takes `--barcodes-fasta` for the bundle but never passes it
   to `render()`. This is the concrete way the "complete CLI surface"
   claim breaks — wire it, and note it in 43a's outcomes.

### 2.2 New params

`params.facility` and `params.analyst_name`, both defaulting to `null`
in `nextflow.config`, threaded into the render context via the params
JSON `collate.nf` already emits (43a §4.4). **Unset renders as `-`**,
never as a blank row that could read as a missing result (rule 18).
Add both to the schema/validation surface `validateParams()` covers if
it enumerates known params.

---

## 3. Tab contents

Per spec 06a-reports.md §6a.2. `render_subjective_input(n)` goes on
every tab, per the boilerplate's existing convention.

| # | Tab | Source in `metadata.json` |
|---|---|---|
| 1 | **Overview** (default) | Key findings (43a's generator), `warnings` mirror, sample metadata + params/versions modals, run context + wall time, `read_qc` |
| 2 | **Assembly** | `assembly.*`, `assembly.bin_metadata`, `annotation`, Bandage PNG, organelle map SVG, annotation GFF |
| 3 | **Validation** | `coverage.recruitment`, `coverage.gate`, `coverage.estimate`, `homology.top_hits[]`, `annotation.cds_crosscheck`, `annotation.genetic_code_agreement` |
| 4 | **Barcodes** | `barcodes.loci[]`, `barcodes.n_passed`, `barcodes.fasta` |

**Every warning raised in a tab must also be appended to 43a's
`warnings` context list**, so it appears in the Overview mirror. A tab
keeps its warning in context — a reader on Validation must see the
coverage caveat there — but no warning may exist *only* behind a tab
(43a §6, CONSTITUTION.md principle 7).

### 3.1 Tab 1 — Overview

Four blocks in order: Key findings, warnings mirror, inputs + run
context, sequencing data quality.

Read QC is a summary table of raw vs clean side by side plus the derived
`filter_yield`. Embed NanoPlot's own `NanoPlot-report.html` as a base64
`data:` URI behind a "full read QC report" button — 43a passes it in.
**Do not** re-render NanoPlot's distribution plots from its standalone
Plotly HTML files; that is fragile against every NanoPlot release
(rule 19), and 43a §4.1 already records the decision.

Read QC sits on Overview rather than an analysis tab because it is an
*input* property, not a result: it tells the reader what they submitted
and how much survived filtering.

### 3.2 Tab 2 — Assembly

Assembly statistics (contigs, N50, total bp) and the per-contig
breakdown (length, coverage, circularised y/n) from
`assembly.contigs[]`. **No separate polishing section** — Flye's
built-in single-iteration pass is the baseline until MEDAKA is
un-deferred (spec §7 Q6). **No SPAdes wording anywhere** (§6a.3).

The Bandage graph image, the per-contig mean coverage chart (§4.1), the
binning classification (target / secondary / off-target) from
`assembly.bin_metadata`, and the plastid QC signals of §5.6.

**No BUSCO.** The mock included it; our spec does not (§6a.3).

The annotated organelle map lives here too, as the visual form of the
same assembly: inline SVG, injected into the DOM rather than via
`<img>` so task 44's per-feature tooltips and click targets work.
Tooltip carries name, product, strand, coordinates, feature source
(`miniprot` vs `mitos`) and the confidence triplet — under §5.4's rules,
which govern how that triplet may be shown.

The annotation GFF downloads from this tab (spec §7).

### 3.3 Tab 3 — Validation

Recruitment + coverage gate (§5.1), homology (§5.2), and the annotation
cross-checks (§5.5), in that order. When a sample soft-failed, **this
tab is the terminal content** — the Assembly and Barcodes tabs must say
no assembly was attempted rather than render empty tables.

### 3.4 Tab 4 — Barcodes

Per-locus recovery table, drop-out reasons (§5.8), and the FASTA
download (§2.1 item 4). Individual sequences are viewable in-page via
the boilerplate's existing `sequence-display.html` macro and modal —
reuse them rather than building new markup.

CDS features come from the broad miniprot pass for every target — the
same features the Assembly tab's map is drawn from, so a locus shown
there and a locus shipped in `barcodes.fasta` are **the same object**.
Say so; a reader should not have to guess whether they are looking at
two views of one thing or two independent analyses.

---

## 4. Charts

### 4.1 "Coverage-along-genome" is not computable — descope it

spec §6a.2's assembly row asks for a per-contig coverage plot. **No
stage in this pipeline computes per-position depth.** Flye's
`assembly_info.txt` gives one mean coverage figure per contig, and that
is what `assembly.contigs[]` carries. Producing a real depth track means
a new alignment + `samtools depth` stage — a pipeline change, not a
reporting one.

**Decision:** render **per-contig mean coverage** as a bar chart, with
target / secondary / off-target binning distinguished by colour. That is
what the section is actually for — the reader needs to see the target
contig sitting at a different depth from the off-target ones — and it is
fully served by data in hand. Spec §6a.2 already says per-contig mean
coverage; do **not** reintroduce a dangling promise of a depth track.
Add a `tasks/todo.md` entry recording that a real depth track is a new
stage and a new task.

### 4.2 Whether Plotly stays

spec §6a.1 says keep Plotly "unless it turns out to bloat the file
substantially at multi-sample scale" — so **measure it, and record the
number in §9** rather than deciding by feel.

`plotly-basic-3.0.0.min.js` is 1.0 MB; with jQuery, Bootstrap and
DataTables the inlined static payload is ~1.6 MB per report, in every
per-sample report *and* the run report. 43a's real integration run
produced a 1.83 MB `report.html` for `INT-ANIMAL-01` with only the
scaffold rendered, so the payload dominates the file today. At that size
a 50-sample batch is ~90 MB of reports — probably acceptable for a
biosecurity archive, but make the call explicitly.

Charts actually needed are all bar charts: filter yield (Overview),
per-contig coverage (Assembly), possibly barcode recovery (Barcodes). If
the payload is judged too heavy, the lever is `EXCLUDE_JS` — the
boilerplate already has the hook — or dropping to inline SVG bars.

**Do not** switch to a CDN-loaded Plotly. §6a.1's self-contained-offline
requirement is not negotiable, and an archived report that renders blank
in three years when the CDN moves is exactly what rule 18 forbids.

---

## 5. What the report must not get wrong

This section is the reason the task exists. Each item is a specific way
a reader is misled, and each is testable.

### 5.1 `low_coverage` is a warned partial result, not a failure

The Validation tab carries a prominent warning that **the assembly was
produced under the warn floor**, and that missing genes and fragmented
contigs are **expected** — the reader must not read incompleteness as a
negative finding (principle 7's second half; spec 02-stages.md §2.1.3).

State the estimate's optimism explicitly (spec §2.1.6): actual depth on
the assembled contig is typically well below `estimated_cov`. On plant
targets, show the sibling-organelle split — `target_assigned_bases`,
`sibling_assigned_bases`, `sibling_organelle_fraction` — so the reader
sees what the estimate was corrected for (§2.1.5). Where `coverage_basis`
is `total_recruited`, **say so**: the sibling split was unavailable and
the estimate may over-state target depth. This is not a footnote —
every real fixture currently reports `total_recruited`.

Where `subsampled` is true, show the seqtk `fraction`, the `seed`, and
`pre_subsample_cov` vs `post_subsample_cov`.

### 5.2 Homology below species threshold

The Validation tab notes explicitly when top-hit identity is below the
species threshold. This **flags** a novel or under-represented taxon; it
does not *assign* one. Taxonomic assignment is Taxodactyl's job
(brief.md §2, principle 1) and no wording in this tab may imply
otherwise. Their proposed-taxonomy panel stays in the mock (§6a.3).

### 5.3 `annotator_failed` renders distinctly from `ok_cds_only`

Both ship a CDS-only annotation. `ok_cds_only` is the honest provisional
state of a target with no annotator configured. `annotator_failed` is a
configured annotator that ran and returned nothing. Show
`annotation.reason`, the tool name and `annotation.annotator_exit_code` —
not just the CDS-only feature set. **A reader must never mistake a tool
failure for a design decision** (task 41).

### 5.4 Confidence is two independent axes, not one score

`annotation.cds_scores[]` carries `pident`, `qcovhsp` and `bitscore` per
feature. These are two axes and must be presented as two:

- **High identity + low coverage** = a truncated fragment of a
  well-referenced gene.
- **Low identity + high coverage** = a complete gene in an
  under-referenced clade — **and this must not read as failure.** This is
  the biosecurity case that matters most, and it is not hypothetical:
  `INT-ANIMAL-01`'s real output has `COX3` at `pident` 39.5 with
  `qcovhsp` 96, and `ATP6` at 39.3 / 100. A report that collapses those
  to a single red "39% confidence" badge is actively misleading about an
  intact gene.
- A **`null` triplet** (no blastp hit against its own gene's panel)
  renders as **explicitly unscored**, not as zero confidence.

Scores may **sort or style** features and nothing else — never which
features appear, never sample status (task 40 §5).

**The above is stated negatively, which is not enough to build against.**
Four positive rules; the first is the load-bearing one:

1. **No composite score exists anywhere** — not in the template, not on
   the result class, not in the rendered DOM. `pident` and `qcovhsp` are
   never averaged, multiplied, weighted, or reduced to a single number,
   badge, bar or colour. This is the invariant that actually closes the
   risk, and unlike "must not read as failure" it is mechanically
   checkable (§6).
2. **Two values, two labels, always adjacent.** Every surface that shows
   one shows the other, separately labelled — the `cds_scores` table
   columns and the map's SVG tooltip alike. If colour is used at all it
   encodes **one** axis; the other uses a different channel (bar length,
   position). Encoding both in one channel is rule 1 by another route.
3. **Name the quadrant affirmatively.** Derive a categorical
   interpretation on the result class (not in the template), with four
   outcomes carrying their own words, e.g. `complete_well_referenced`,
   `truncated_well_referenced`, `complete_under_referenced`,
   `fragmentary`. The low-identity/high-coverage case gets an
   **affirmative** label — "complete gene, under-referenced clade" — not
   merely the absence of a warning. A reader must be *told* this is a
   good result, not left to infer it from the lack of a red badge.
   Thresholds separating the quadrants are constants at the top of the
   module, not magic numbers in a template.
4. **Show the interpretation key in the tab**, in the two-axis language
   of the first two bullets, so the reader learns the axes rather than
   guessing at them.

**`bitscore` must not sort features across different genes.** It scales
with alignment length, so ranking a feature set by it puts long genes
above short ones irrespective of quality — a `COX1` fragment outranking
an intact `ATP8`. Use it only to compare hits for the *same* gene, or not
at all. `INT-ANIMAL-01`'s real spread makes this concrete: `COX3` 196.0
vs `ND3` 66.2, where the `ND3` call is the better-covered one at 85%
`qcovhsp`.

### 5.5 Surface the QC disagreements the merge would otherwise absorb

Two signals from `annotation`, on the Validation tab (tasks/todo.md,
task 31 carry-forward). **43a already raises both as `warnings` entries**
— this task renders the detail the warning points at, and must not
duplicate or contradict the warning text:

- **`cds_crosscheck`** — `miniprot_only`, `annotator_only` and
  `coordinate_conflicts`. Two independent methods disagreeing on a gene
  is a real QC signal. `INT-ANIMAL-01` has two `annotator_only` entries.
  Note the gene-name convention task 42 §5.3 settled: these carry both
  `raw` (the annotator's fragment-suffixed strings, `atp8_1`) and
  `canonical` (`ATP8`) lists. Display canonical, expose raw.
- **`genetic_code_agreement: false`** — ANNOTATE and EXTRACT_BARCODES
  trialled different tables. Compare `genetic_code_annotate` against
  `genetic_code_cds` and say so when they differ.

Delete the task 31 todo.md entry once done.

### 5.6 Withheld plastid substitution

tasks/todo.md, task 24 carry-forward. The Assembly tab renders, from
`assembly.bin_metadata.plastid_canonicalisation`:

- the **edge-count reading** (`edge_count`, `lsc_edge`, `ir_edge`,
  `ssc_edge`);
- the **"canonical structure, no taxonomic support"** warning when
  `substitution_applied` is false, with `non_canonical_reason` /
  `substitution_withheld_reason`;
- the **`target_source`** provenance line when it was applied
  (spec 03-organelles.md §3.6).

`INT-PLANT-01-pt` exercises the applied branch (`substitution_applied:
true`, `edge_count: 3`). The withheld branch has no real fixture — cover
it in unit tests (§6).

Delete the task 24 todo.md entry once done.

### 5.7 Plant arms are CDS-only and must say so

tRNA/rRNA come from MITOS2 on `animal_mt` only. On `plant_pt` and
`plant_mt` the Assembly and Barcodes tabs **state** that the annotation
is CDS-only pending spec §8 item 3, rather than leaving the absence of
tRNA/rRNA to be inferred from an empty row (principle 7). Both plant
fixtures currently report `ok_cds_only`, so this renders on two of the
three real samples.

### 5.8 Barcode drop-outs carry their reason

The Barcodes tab lists panel loci **not** recovered, each with its
reason — `not_found`, `invalid_length`, `identity_below_floor`,
`internal_stop_codon` — from `barcodes.loci[]`. An absent locus rendered
as a blank row is the "output indistinguishable from a confident
negative" principle 7 forbids. `INT-ANIMAL-01` is a live example: `COX1`
fails on `internal_stop_codon` while the sample overall is `ok`.

---

## 6. Tests — extending `scripts/tests/test_report.py`

Rule 14 applies: **100% branch coverage** via `scripts/pytest.sh`. Mock
at the tool boundary only (task 27) — build real `metadata.json`
fixtures, render real HTML, assert on the real output string.

Scaffold corrections (§2.1):

- **The H1 title tracks `organelle`** — `mt` renders "Mitochondrial",
  `pt` renders "Chloroplast", and neither renders the other's word.
- **Key findings is inside the Overview pane**, not in the
  always-visible header — assert on DOM position, since a test that only
  checks the string is present passes either way.
- **Overview is the default-active tab.** This is the principle-7
  guarantee now that Key findings is no longer always visible; if a
  future edit makes another tab active on load, negative clarity is
  behind a click and the test must fail.
- **Wall time renders**, and **unset facility/analyst render `-`**, not
  a blank cell (§2.2).
- **The barcode FASTA download is present** and carries the sequences
  from `barcodes.fasta`.

Content rules (§5):

- **`low_coverage`** renders the warn-floor warning in **both** the
  Validation tab and the Overview mirror (§5.1) — the mirror assertion is
  what stops the tab layout from swallowing a principle-7 signal.
- **`coverage_basis: "total_recruited"`** renders the over-statement
  caveat.
- **`subsampled: true`** renders fraction, seed and pre/post coverage.
- **A soft-failed sample** renders Validation as terminal, with Assembly
  and Barcodes stating no assembly was attempted rather than showing
  empty tables.
- **Top-hit identity below the species threshold** renders the flag, and
  no wording asserting a taxonomic assignment (§5.2).
- **`annotator_failed`** renders `reason` and `annotator_exit_code`, and
  produces different text from `ok_cds_only` (§5.3).
- **A `null` score triplet** renders as unscored, and the string "0" does
  not appear as its confidence (§5.4).
- **No composite score exists** (§5.4 rule 1) — the regression that
  matters most, and the one assertion here that cannot be satisfied
  vacuously. Assert the result class exposes no combined/mean/weighted
  score attribute, and that a feature's rendered row contains `pident`
  and `qcovhsp` as two separately-labelled values.
- **The quadrant classifier is unit-tested directly**, at and either side
  of each threshold, independently of any rendering (§5.4 rule 3).
- **A low-identity / high-coverage feature carries its affirmative
  label** — assert the "under-referenced" wording is *present*, not
  merely that a danger class is absent. Use `COX3` at 39.5 / 96 from
  `INT-ANIMAL-01` as the fixture, so the test fails if a future edit
  reintroduces a single-scale reading of a real intact gene.
- **`bitscore` does not drive cross-gene sort order** (§5.4) — assert a
  feature set containing a high-bitscore long gene and a better-covered
  short gene does not rank purely by bitscore.
- **`cds_crosscheck` non-empty** and **`genetic_code_agreement: false`**
  both surface (§5.5).
- **`substitution_applied: false`** renders the "canonical structure, no
  taxonomic support" warning; **`true`** renders the `target_source`
  line (§5.6).
- **`ok_cds_only` on a plant target** states CDS-only explicitly (§5.7).
- **Every barcode drop-out reason** renders its reason string (§5.8).
- **Self-containment holds** after every tab is added — re-assert 43a's
  no-external-fetch check, since this task adds the most markup. Note
  43a's two allowlisted strings (the repo subtitle link and the inert
  `plotly.com` attribution baked into the vendored JS); do not widen the
  allowlist further without recording why.

---

## 7. Integration reconciliation

- **`tests/integration/assertions.sh`.** Extend 43a's block: for
  `INT-PLANT-01-mt` (the only real `low_coverage` fixture, at 28.48×)
  assert the warn-floor warning string is present in `report.html`; for
  `INT-PLANT-01-pt` assert the plastid canonicalisation block rendered;
  for `INT-ANIMAL-01` assert the `COX1` `internal_stop_codon` drop-out
  reason appears.
- **Assert on structure and status, not on numbers** that shift with
  reference-data rebuilds (rule 19). Assert the coverage figure is
  *present and labelled*, never that it equals 28.48.
- **Stub run stays green**, including the new `workflow.start` argument
  and the `barcodes.fasta` renderer input.
- **Known gap:** no fixture exercises `fail`, `no_assembly` or
  `no_barcode` with real data (task 42 §9). §6's unit tests remain the
  only coverage of the soft-fail tab behaviour. Record it; do not paper
  over it.

---

## 8. Spec updates

spec/06a-reports.md §6a.1, §6a.2 and §6a.5 were rewritten to the
four-tab structure ahead of this task (2026-09-03), so the tab taxonomy,
the organelle-specific title, the wall-time params and the barcode FASTA
download are already specified. Remaining reconciliation:

- **tasks/completed/43a_report_scaffold.md** — append to §14 Outcomes
  that criteria 5 and 7 are superseded by this task (§2), with the
  reason. Do not silently leave a completed task claiming a contract
  this task changes.
- **tasks/todo.md** — delete the task 24 per-sample-report entry (§5.6)
  and the task 31 `cds_crosscheck` entry (§5.5), both resolved here. Add
  the depth-track entry from §4.1. Leave the `ORGANELLE_MAP` and
  `RUN_REPORT` entries alone.

Per the task-authoring convention: when editing spec files, reference
tasks as plain text (`task 43b_report_stage_tabs.md`), never as markdown
links, and never with line anchors.

---

## 9. Acceptance criteria

1. The strip carries exactly four tabs — Overview, Assembly, Validation,
   Barcodes — with Overview active on load.
2. All four panes render from `metadata.json` plus the file-shaped
   renderer arguments; flake8-clean at 79 columns with 100% branch
   coverage under `scripts/pytest.sh`.
3. The four §2.1 scaffold corrections are in place: organelle-specific
   H1, Key findings inside Overview, wall-time component restored with
   its two new params, `barcodes.fasta` wired to the renderer.
4. Every §5 item is implemented and individually covered by a test, with
   §5.4 rule 1 (no composite score) asserted mechanically.
5. Every tab warning also appears in the Overview mirror (§3).
6. The report remains a **single self-contained file** with no external
   fetch after all tabs are added.
7. The annotation GFF and the barcode FASTA both download from the
   report as embedded `data:` URIs.
8. A soft-failed sample renders Validation as terminal content, not as
   empty downstream tables (§5.1).
9. Plotly's size impact measured and the keep/drop decision recorded
   (§4.2).
10. `-stub-run` green; `assertions.sh` extended and passing on the
    nightly integration profile.
11. Spec, 43a outcomes and todo.md reconciled per §8.

---

## 10. Out of scope

- **The real organelle map SVG** — task 44_organelle_map.md. Inline
  whatever is given, including the current zero-byte stub.
- **`run-report.html` / `run_manifest.json`** — task 45_run_report.md.
- **A per-position depth track** — §4.1; a new stage, a new task.
- **Adding `INT-PLANT-01-mt` to `ASSEMBLING_SAMPLES`** — task
  36_plant_mt_integration_coverage.md, shelved.
- **Plant annotation arms** — P4b. This task *reports* `ok_cds_only`
  (§5.7); it does not improve it.

---

## 11. Outcomes

<!-- Filled in on completion. -->
