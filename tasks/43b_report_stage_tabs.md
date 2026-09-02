# Task 43b — Per-sample `report.html`: the seven stage tabs

**Phase:** P4a (spec 06-phases.md).
**Goal:** Fill in tabs 2–8 of the per-sample report — one per pipeline
stage group, in pipeline order — on the scaffold task
43a_report_scaffold.md ships.

**Depends on:** task 43a_report_scaffold.md. It fixes the renderer's
complete CLI surface, the eight-tab strip, the `warnings` mirror
mechanism, the always-visible Key findings header, and every plumbing
input including `metadata.json`'s `read_qc` block. **Do not start before
it lands** — this task should be filling panes, not restructuring the
frame.

**Related tasks:** 44_organelle_map.md supplies the real SVG for tab 8;
this task builds against its zero-byte stub and inlines whatever it is
given, so task 44's tool decision is not a blocker. 45_run_report.md
reuses the shared components; factor rather than inline.

---

## 1. Overview

43a answered "did this sample work". This task answers "how confident
should I be" and "what exactly do I have".

Almost all of it is presentation of data already in `metadata.json`.
Very little new logic is needed. What makes this task worth a brief of
its own is §4: seven specific ways a reader can be misled by a correct
number rendered carelessly, each traceable to a CONSTITUTION principle
and each individually testable. That section is the task; the tab
layout is the easy part.

---

## 2. Tab contents

Section order follows spec 06a-reports.md §6a.2, which is pipeline order.
Tab 1 (Overview) is 43a's.

| # | Tab | Source in `metadata.json` |
|---|---|---|
| 2 | Sequencing data quality | `read_qc` (43a §4.1) |
| 3 | Recruitment + coverage gate | `coverage.gate`, `coverage.estimate`, recruit stats |
| 4 | *De novo* assembly | `assembly.contig_count` / `n50` / `total_bp` / `contigs[]` |
| 5 | Assembly quality | `assembly.bin_metadata`, Bandage PNG |
| 6 | Homology | `homology.top_hits[]` |
| 7 | Barcode panel | `barcodes.loci[]`, `barcodes.n_passed` |
| 8 | Annotated organelle map | `annotation`, organelle map SVG |

`render_subjective_input(n)` goes on every tab, per the boilerplate's
existing convention.

**Every warning raised in a stage tab must also be appended to 43a's
`warnings` context list**, so it appears in the Overview mirror. A stage
tab keeps its warning in context — a reader on tab 3 must see the
coverage caveat there — but no warning may exist *only* behind a tab
(43a §6, CONSTITUTION.md principle 7).

### 2.1 Tab 2 — read QC

Summary table of raw vs clean side by side, plus the derived
`filter_yield`. Embed NanoPlot's own `NanoPlot-report.html` as a base64
`data:` URI behind a "full read QC report" button — 43a passes it in.
**Do not** re-render NanoPlot's distribution plots from its standalone
Plotly HTML files; that is fragile against every NanoPlot release
(rule 19), and 43a §4.1 already records the decision.

### 2.2 Tab 4 — assembly

Assembly statistics (contigs, N50, total bp) and the per-contig
breakdown (length, coverage, circularised y/n) from
`assembly.contigs[]`. **No separate polishing section** — Flye's
built-in single-iteration pass is the baseline until MEDAKA is
un-deferred (spec §7 Q6). **No SPAdes wording anywhere** (§6a.3).

### 2.3 Tab 5 — assembly quality

The Bandage graph image, the per-contig coverage chart (§3.1), the
binning classification (target / secondary / off-target) from
`assembly.bin_metadata`, and the plastid QC signals of §4.6.

**No BUSCO.** The mock included it; our spec does not (§6a.3).

### 2.4 Tab 8 — annotated organelle map

Inline SVG, injected into the DOM rather than via `<img>` so task 44's
per-feature tooltips and click targets work. Tooltip carries name,
product, strand, coordinates, feature source (`miniprot` vs `mitos`) and
the confidence triplet — under §4.4's rules, which govern how that
triplet may be shown.

CDS features come from the broad miniprot pass for every target — the
same features the barcode panel in tab 7 is drawn from, so a locus shown
here and a locus shipped in `barcodes.fasta` are **the same object**.
Say so; a reader should not have to guess whether they are looking at two
views of one thing or two independent analyses.

---

## 3. Charts

### 3.1 "Coverage-along-genome" is not computable — descope it

spec §6a.2's assembly-quality row asks for a "per-contig
coverage-along-genome plot". **No stage in this pipeline computes
per-position depth.** Flye's `assembly_info.txt` gives one mean coverage
figure per contig, and that is what `assembly.contigs[]` carries.
Producing a real depth track means a new alignment + `samtools depth`
stage — a pipeline change, not a reporting one.

**Decision:** render **per-contig mean coverage** as a bar chart, with
target / secondary / off-target binning distinguished by colour. That is
what the section is actually for — the reader needs to see the target
contig sitting at a different depth from the off-target ones — and it is
fully served by data in hand. Amend spec §6a.2 to say per-contig
coverage, and do **not** leave a dangling promise of a depth track. Add a
`tasks/todo.md` entry recording that a real depth track is a new stage
and a new task.

### 3.2 Whether Plotly stays

spec §6a.1 says keep Plotly "unless it turns out to bloat the file
substantially at multi-sample scale" — so **measure it, and record the
number in §8** rather than deciding by feel.

`plotly-basic-3.0.0.min.js` is 1.0 MB; with jQuery, Bootstrap and
DataTables the inlined static payload is ~1.6 MB per report, in every
per-sample report *and* the run report. At that size a 50-sample batch is
~80 MB of reports — probably acceptable for a biosecurity archive, but
make the call explicitly.

Charts actually needed are all bar charts: filter yield (tab 2),
per-contig coverage (tab 5), possibly barcode recovery (tab 7). If the
payload is judged too heavy, the lever is `EXCLUDE_JS` — the boilerplate
already has the hook — or dropping to inline SVG bars.

**Do not** switch to a CDN-loaded Plotly. §6a.1's self-contained-offline
requirement is not negotiable, and an archived report that renders blank
in three years when the CDN moves is exactly what rule 18 forbids.

---

## 4. What the report must not get wrong

This section is the reason the task exists. Each item is a specific way a
reader is misled, and each is testable.

### 4.1 `low_coverage` is a warned partial result, not a failure

Tab 3 carries a prominent warning that **everything below it was
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

When a sample soft-failed, **tab 3 is the terminal content** — no
assembly was attempted, and tabs 4–8 must say that rather than render
empty tables.

### 4.2 Homology below species threshold

Tab 6 notes explicitly when top-hit identity is below the species
threshold. This **flags** a novel or under-represented taxon; it does not
*assign* one. Taxonomic assignment is Taxodactyl's job (brief.md §2,
principle 1) and no wording in this tab may imply otherwise. Their
proposed-taxonomy panel stays in the mock (§6a.3).

### 4.3 `annotator_failed` renders distinctly from `ok_cds_only`

Both ship a CDS-only annotation. `ok_cds_only` is the honest provisional
state of a target with no annotator configured. `annotator_failed` is a
configured annotator that ran and returned nothing. Show
`annotation.reason`, the tool name and `annotation.annotator_exit_code` —
not just the CDS-only feature set. **A reader must never mistake a tool
failure for a design decision** (task 41).

### 4.4 Confidence is two independent axes, not one score

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
   checkable (§5).
2. **Two values, two labels, always adjacent.** Every surface that shows
   one shows the other, separately labelled — the `cds_scores` table
   columns and tab 8's SVG tooltip alike. If colour is used at all it
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

### 4.5 Surface the QC disagreements the merge would otherwise absorb

Two signals from `annotation` (tasks/todo.md, task 31 carry-forward):

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

### 4.6 Withheld plastid substitution

tasks/todo.md, task 24 carry-forward. Tab 5 renders, from
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
it in unit tests (§5).

Delete the task 24 todo.md entry once done.

### 4.7 Plant arms are CDS-only and must say so

tRNA/rRNA come from MITOS2 on `animal_mt` only. On `plant_pt` and
`plant_mt` tabs 7 and 8 **state** that the annotation is CDS-only pending
spec §8 item 3, rather than leaving the absence of tRNA/rRNA to be
inferred from an empty row (principle 7). Both plant fixtures currently
report `ok_cds_only`, so this renders on two of the three real samples.

### 4.8 Barcode drop-outs carry their reason

Tab 7 lists panel loci **not** recovered, each with its reason —
`not_found`, `invalid_length`, `identity_below_floor`,
`internal_stop_codon` — from `barcodes.loci[]`. An absent locus rendered
as a blank row is the "output indistinguishable from a confident
negative" principle 7 forbids. `INT-ANIMAL-01` is a live example: `COX1`
fails on `internal_stop_codon` while the sample overall is `ok`.

---

## 5. Tests — extending `scripts/tests/test_report.py`

Rule 14 applies: **100% branch coverage** via `scripts/pytest.sh`. Mock
at the tool boundary only (task 27) — build real `metadata.json`
fixtures, render real HTML, assert on the real output string.

- **`low_coverage`** renders the warn-floor warning in **both** tab 3 and
  the Overview mirror (§4.1) — the mirror assertion is what stops the tab
  layout from swallowing a principle-7 signal.
- **`coverage_basis: "total_recruited"`** renders the over-statement
  caveat.
- **`subsampled: true`** renders fraction, seed and pre/post coverage.
- **A soft-failed sample** renders tab 3 as terminal, with tabs 4–8
  stating no assembly was attempted rather than showing empty tables.
- **Top-hit identity below the species threshold** renders the flag, and
  no wording asserting a taxonomic assignment (§4.2).
- **`annotator_failed`** renders `reason` and `annotator_exit_code`, and
  produces different text from `ok_cds_only` (§4.3).
- **A `null` score triplet** renders as unscored, and the string "0" does
  not appear as its confidence (§4.4).
- **No composite score exists** (§4.4 rule 1) — the regression that
  matters most, and the one assertion here that cannot be satisfied
  vacuously. Assert the result class exposes no combined/mean/weighted
  score attribute, and that a feature's rendered row contains `pident`
  and `qcovhsp` as two separately-labelled values.
- **The quadrant classifier is unit-tested directly**, at and either side
  of each threshold, independently of any rendering (§4.4 rule 3).
- **A low-identity / high-coverage feature carries its affirmative
  label** — assert the "under-referenced" wording is *present*, not
  merely that a danger class is absent. Use `COX3` at 39.5 / 96 from
  `INT-ANIMAL-01` as the fixture, so the test fails if a future edit
  reintroduces a single-scale reading of a real intact gene.
- **`bitscore` does not drive cross-gene sort order** (§4.4) — assert a
  feature set containing a high-bitscore long gene and a better-covered
  short gene does not rank purely by bitscore.
- **`cds_crosscheck` non-empty** and **`genetic_code_agreement: false`**
  both surface (§4.5).
- **`substitution_applied: false`** renders the "canonical structure, no
  taxonomic support" warning; **`true`** renders the `target_source`
  line (§4.6).
- **`ok_cds_only` on a plant target** states CDS-only explicitly (§4.7).
- **Every barcode drop-out reason** renders its reason string (§4.8).
- **Self-containment holds** after every tab is added — re-assert 43a's
  no-external-fetch check, since this task adds the most markup.

---

## 6. Integration reconciliation

- **`tests/integration/assertions.sh`.** Extend 43a's block: for
  `INT-PLANT-01-mt` (the only real `low_coverage` fixture, at 28.48×)
  assert the warn-floor warning string is present in `report.html`; for
  `INT-PLANT-01-pt` assert the plastid canonicalisation block rendered;
  for `INT-ANIMAL-01` assert the `COX1` `internal_stop_codon` drop-out
  reason appears.
- **Assert on structure and status, not on numbers** that shift with
  reference-data rebuilds (rule 19). Assert the coverage figure is
  *present and labelled*, never that it equals 28.48.
- **Stub run stays green.**
- **Known gap:** no fixture exercises `fail`, `no_assembly` or
  `no_barcode` with real data (task 42 §9). §5's unit tests remain the
  only coverage of the soft-fail tab behaviour. Record it; do not paper
  over it.

---

## 7. Spec updates

- **spec/06a-reports.md §6a.2** — change "coverage-along-genome" to
  per-contig mean coverage (§3.1).
- **tasks/todo.md** — delete the task 24 per-sample-report entry (§4.6)
  and the task 31 `cds_crosscheck` entry (§4.5), both resolved here. Add
  the depth-track entry from §3.1. Leave the `ORGANELLE_MAP` and
  `RUN_REPORT` entries alone.

Per the task-authoring convention: when editing spec files, reference
tasks as plain text (`task 43b_report_stage_tabs.md`), never as markdown
links, and never with line anchors.

---

## 8. Acceptance criteria

1. All seven stage tabs render from `metadata.json` plus 43a's
   file-shaped arguments; flake8-clean at 79 columns with 100% branch
   coverage under `scripts/pytest.sh`.
2. Every §4 item is implemented and individually covered by a test, with
   §4.4 rule 1 (no composite score) asserted mechanically.
3. Every stage-tab warning also appears in the Overview mirror (§2).
4. The report remains a **single self-contained file** with no external
   fetch after all tabs are added.
5. The annotation GFF downloads from the report as an embedded `data:`
   URI (43a §4.3 wires it; this task renders the control).
6. A soft-failed sample renders tab 3 as terminal content, not as empty
   downstream tables (§4.1).
7. Plotly's size impact measured and the keep/drop decision recorded
   (§3.2).
8. `-stub-run` green; `assertions.sh` extended and passing on the nightly
   integration profile.
9. Spec and todo.md reconciled per §7.

---

## 9. Out of scope

- **The scaffold, plumbing and Key findings block** — task
  43a_report_scaffold.md. If an input turns out to be unreachable, that
  is a 43a defect; fix it there rather than working around it here.
- **The real organelle map SVG** — task 44_organelle_map.md. Inline
  whatever is given, including the current zero-byte stub.
- **`run-report.html` / `run_manifest.json`** — task 45_run_report.md.
- **A per-position depth track** — §3.1; a new stage, a new task.
- **Adding `INT-PLANT-01-mt` to `ASSEMBLING_SAMPLES`** — task
  36_plant_mt_integration_coverage.md, shelved.
- **Plant annotation arms** — P4b. This task *reports* `ok_cds_only`
  (§4.7); it does not improve it.

---

## 10. Outcomes

<!-- Filled in on completion. -->
