# Task 43 — Stage 15 `COLLATE`: per-sample `report.html`

> **STUB — expand before executing.**

**Phase:** P4a (from spec §6).
**Depends on:** task 42 — the
`metadata.json` contract this renders from.
**Goal:** Wire `scripts/report/boilerplate/`
into a real report image and render the per-sample `report.html` from
`COLLATE`'s bundle.

The boilerplate is already in the repo and entirely unused: ~550 lines of
Python (`report.py`, `config.py`, `results.py`) plus vendored Bootstrap 5,
Plotly, DataTables and the DAFF logo. It ships stubbed by design —
concrete result classes, config accessors and stage-per-tab templates are
added as each stage lands
(spec §6a.5).
This is the task that lands nearly all of them at once, so expect it to be
the bulk of P4a's remaining work.

## Scope sketch

- Section-by-section build-out against
  spec §6a.2's
  table — one H2 per pipeline stage, in pipeline order.
- Adopt from the mock: two-column header (*Key findings* + *Input
  parameters*), modal detail views, interactive charts, inline SVG gene
  map (spec §6a.1).
- **Leave in the mock**
  (spec §6a.3):
  SPAdes stats, the BUSCO section, `core_nt`, their proposed-taxonomy
  panel, their fixed gene panel. The mock informs *presentation only* —
  its scientific spec must not leak in.
- Self-contained single-file HTML, viewable offline, no external fetch.

## Things the report must not get wrong

- **Negative clarity front-and-centre** (brief.md §3.8,
  principle 7). No assembly / no extractable locus
  is stated in *Key findings* at the top, not buried in a later section.
- **`low_coverage` reads as a warned partial result, not a failure.**
  The recruitment/coverage section carries a prominent warning that
  everything below it was produced under the warn floor, and that missing
  genes and fragmented contigs are **expected** here. State the
  estimate's optimism explicitly
  (spec §2.1.6).
- **`annotator_failed` renders distinctly from `ok_cds_only`** (task 41).
  Both ship a CDS-only annotation, but one is a provisional-by-design
  state and the other is a configured tool that ran and returned nothing —
  show the `reason`, tool name and exit code. A reader must never mistake
  a tool failure for a design decision.
- **Confidence is two independent axes, not one score** (task 40).
  High identity + low coverage = truncated fragment of a well-referenced
  gene; low identity + high coverage = complete gene in an
  under-referenced clade — and **the second must not read as failure**
  (this is the biosecurity case that matters most). A `null` triplet
  renders as *explicitly unscored*, not as zero confidence. Scores may
  sort or style features and nothing else — never which features appear,
  never sample status (task 40 §5).
- **Surface `cds_crosscheck` disagreements and
  `genetic_code_agreement: false`** (todo.md, task 31 carry-forward) —
  two independent methods disagreeing on a gene, or ANNOTATE and
  EXTRACT_BARCODES trialling different tables, are real QC signals that
  must not be silently absorbed by the merge.
- **Plant arms are CDS-only and must say so** rather than leaving the
  absence of tRNA/rRNA to be inferred.
- **Withheld plastid substitution** (todo.md, task 24 carry-forward): the
  "Assembly quality assessment" section renders the "canonical structure,
  no taxonomic support" warning and the `target_source` provenance line
  (spec §3.6).
- Annotation GFF available for download (spec §6a "Specific feature
  requests").

## Open questions for the expansion pass

- Does Plotly bloat the file unacceptably at multi-sample scale? Spec
  says keep it unless it does — measure.
- The gene-map SVG arrives from task 44; render
  against the stub SVG until it lands.
