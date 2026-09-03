## 6a. Report design (per-sample `report.html`)

**Design reference: `bin/report/` (task 43a_report_scaffold.md §2).**
`mock-report.html` has never existed in this repository, and every
`../wf-report-boilerplate/` path below is stale — that module moved to
`bin/report/` (Python package) + `scripts/report/templates/` +
`scripts/report/static/` (task 43a §5.1). The boilerplate module *is*
the design reference: it is trimmed from a prior Nextflow workflow
specifically to preserve the paradigms this project extends, so wf5's
report stays structured like the other workflows in the family.
§6a.3's list still applies as a **negative constraint** — those
elements must not appear — even without a mock to point at.

### 6a.1 Design elements to adopt

- **Bootstrap-based single-file HTML** — self-contained, no external asset fetch, viewable offline.
- **Organelle-specific H1 title**, derived at render time rather than hard-coded: `organelle: "mt"` → *Mitochondrial genome assembly*, `organelle: "pt"` → *Chloroplast genome assembly*. A report that calls itself "Organelle assembly + barcode recovery" makes the reader work out which organelle they are looking at from the body; the title should tell them. The mapping lives in the renderer as a constant, not in a template.
- **Everything below the heading block lives in a tab.** The heading block — H1 title, sample ID, DAFF logo, subtitle — is the only always-visible content. *Key findings* is the top of the **Overview** tab, not a separate header row above the strip.
- **Four reader-perspective tabs, not one tab per pipeline stage:** *Overview*, *Assembly*, *Validation*, *Barcodes*. The pipeline's stage boundaries are an implementation detail; a biosecurity officer asks "did it work, what did I get, should I trust it, which barcodes can I use" — the tabs answer those questions in that order. Several stages therefore share a tab, and one stage's output (the coverage gate) is deliberately filed under the question it answers rather than the stage that produced it.
- **Overview is tab 1 and opens by default**, so the first thing a reader sees is *Key findings*. It also mirrors every warning raised on any other tab, via a `warnings` list in the render context — so a caveat is never visible only behind a tab a reader might not open.
- **Modal dialogs for detail views** — e.g. full parameter dump, tool versions with commit/build hashes. Keeps the main flow terse; detail one click away.
- **Interactive charts** where a plot beats a table — filter yield (Overview), per-contig mean coverage (Assembly), barcode recovery (Barcodes). All bar charts. Plotly is the mock's choice; keep unless it turns out to bloat the file substantially at multi-sample scale. Read-length and Q-score distributions are **not** re-rendered from NanoPlot's standalone Plotly HTML — scraping trace JSON out of them is fragile against every NanoPlot release ([rule 19](../CONSTITUTION.md)); link out to NanoPlot's own embedded report instead.
- **Inline SVG for the annotated organelle map** — clickable/hoverable gene features with tooltip metadata (name, strand, coordinates, source, confidence score).
- **"Negative clarity" front-and-centre** ([brief.md §3.8](brief.md)) — when no organelle is assembled or no locus is extractable, the *Key findings* block says so explicitly at the top of the default tab, not buried behind a tab the reader must choose to open. `no_assembly` and `no_barcode` must never read as the same outcome, and `low_coverage` must never read as a failure (task 43a §7).

### 6a.2 Tab outline (our pipeline's content)

The heading block (H1 organelle-specific title per §6a.1, sample ID,
logo, subtitle) sits above the strip. Everything else is in one of four
tabs.

| Tab | Content source | Notes |
|---|---|---|
| **1. Overview** (default) | Derived summary + `samples.csv` (§0) + `NANOPLOT_RAW`/`NANOPLOT_CLEAN` + run provenance | Four blocks, in this order. **(a) Key findings** — bulleted plain-English summary: assembly outcome, coverage, top BLAST hit, gene-feature count, panel-locus extraction count. This is the negative-clarity surface of §6a.1's last bullet and it is the first thing on the first tab. **(b) The warnings mirror** — every warning raised on any other tab, each labelled with the tab that owns it. **(c) Inputs and run context** — sample_id, kingdom, organelle, and any submitter-supplied optional columns from `samples.csv`; the "view all parameters" and "view tool versions" modals; the facility / analyst / analysis-started / analysis-completed / wall-time panel; and the run provenance lines (pipeline commit, reference-bundle version and build date). **(d) Sequencing data quality** — pre/post read stats and the filter-yield summary from CHOPPER + FILTLONG, from `metadata.json`'s `read_qc`. Read QC is an *input* property, not a result, which is why it sits here rather than with the analysis tabs: it tells the reader what they submitted and how much of it survived filtering. Link out to NanoPlot's own full report rather than re-rendering its distribution plots. |
| **2. Assembly** | `METAFLYE` + `BANDAGE_NG` + `BIN_TARGET` + `ANNOTATE` → `JOIN_ANNOTATION_SCORES` | "What did I get." Assembly statistics (contigs, N50, total bp) and the per-contig breakdown (length, coverage, circularised y/n). Polishing is Flye's built-in single-iteration pass — no separate polish section until MEDAKA is un-deferred (§7 Q6). Assembly graph image, per-contig mean coverage chart, binning classification (target / secondary / off-target). Carries the plastid QC signals of [§3.6](03-organelles.md#36-plastid-quadripartite-canonicalisation): the edge-count reading, the "canonical structure, no taxonomic support" warning when the `path1` substitution was withheld, and `target_source` provenance when it was applied. **No BUSCO** — mock includes it, our spec does not (kingdom-agnostic organelle BUSCO sets are patchy; revisit if needed). **The annotated organelle map lives here too**, as the visual form of the same assembly: inline SVG, clickable/hoverable gene features. tRNA/rRNA features come from MITOS2 on `animal_mt`; plant arms are CDS-only pending [§8 item 3](07-open-questions.md#8-remaining-open-questions) and **must say so** rather than leaving the absence to be inferred ([principle 7](../CONSTITUTION.md)). Tooltip shows name, product, strand, coordinates, feature source, and confidence: `pident`/`qcovhsp`/`bitscore` from `annotation_summary.json`'s `cds_scores` (task 40), read as two independent axes — high identity + low coverage is a truncated fragment of a well-referenced gene, low identity + high coverage is a complete gene in an under-referenced clade, and the second must not read as failure. A feature with a `null` triplet (no blastp hit against its own gene's panel) renders as explicitly unscored, not as zero confidence. The report may sort or style features by score — that is the only consequence a score is permitted to have (task 40 §5); it never changes which features appear or a sample's status. `status: "annotator_failed"` (task 41) must render distinctly from `ok_cds_only`: both ship a CDS-only annotation, but `ok_cds_only` is the honest provisional state of a target with no annotator configured, while `annotator_failed` is a configured annotator that ran and returned nothing — the `reason` field names the tool and exit code and must be shown, not just the CDS-only feature set. A reader must never mistake a tool failure for a provisional-by-design state. The annotation GFF downloads from here (§7). |
| **3. Validation** | `RECRUIT` + `COVERAGE_GATE` + `BLAST_VALIDATE` + `annotation.cds_crosscheck` | "Should I trust it." Everything bearing on the result's trustworthiness, gathered in one place regardless of which stage produced it. **Recruitment + coverage gate:** recruited-read count and total bases, estimated coverage vs the configured HARD MIN / WARN / MAX (§2.1), gate decision (pass / pass-with-warning / subsampled / soft-fail) with the seqtk fraction and seed if applicable. On a **pass with warning** (`low_coverage`, [§2.1.3](02-stages.md#213-decision-matrix)) this tab carries a prominent warning that the assembly was produced under the warn floor and that missing genes and fragmented contigs are expected — the reader must not read incompleteness as a negative finding. State the estimate's optimism explicitly ([§2.1.6](02-stages.md#216-provenance-of-the-two-floors)): actual depth on the assembled contig is typically well below `estimated_cov`. On plant targets, also the sibling-organelle split of the recruited pool from `sample_status.json` — `target_assigned_bases` vs `sibling_assigned_bases` and the resulting `sibling_organelle_fraction` ([§2.1.5](02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)). Where `coverage_basis` is `total_recruited`, say so: the sibling split was unavailable and the estimate may over-state target depth. **Homology:** top hits per target contig against the kingdom-appropriate RefSeq organelle DB (§4.2) — identity, coverage, subject accession. Explicit note when top-hit identity is below species-threshold — flags novel/underrepresented taxa without pretending to *assign* them (assignment is Taxodactyl's job — see [brief.md §2](brief.md) boundary). **Annotation cross-checks:** where `annotation_summary.json` records `cds_crosscheck` disagreements between miniprot and MITOS2, surface them — two independent methods disagreeing on a gene is a real QC signal — along with `genetic_code_agreement: false`. When a sample soft-failed, **this tab is the terminal content** — no assembly was attempted, and the Assembly and Barcodes tabs must say so rather than render empty tables. |
| **4. Barcodes** | `EXTRACT_BARCODES` + `barcodes.fasta` | "Which barcodes can I use." Per-locus: extracted length, ORF status, genetic code used, coordinates, protein-identity to the miniprot reference. Panel loci **not** recovered are listed with the reason (`not_found` / `invalid_length` / `identity_below_floor` / `internal_stop_codon`) — an absent locus rendered as a blank row is the "output indistinguishable from a confident negative" [principle 7](../CONSTITUTION.md) forbids. **The recovered barcodes download as FASTA from this tab**, embedded as a `data:` URI so the download survives the report being emailed on its own; individual sequences are viewable in-page via the boilerplate's sequence-display modal. The loci here are drawn from the same miniprot CDS features the Assembly tab's map renders, so a locus shown on the map and a locus shipped in `barcodes.fasta` are **the same object** — say so, rather than leaving a reader to guess whether they are two views of one thing or two independent analyses. |

### 6a.3 Elements to leave in the mock

- **SPAdes references / SPAdes-specific stats.** We use metaFlye — swap all assembler mentions.
- **BUSCO section.** Not in our spec. Do not add.
- **`core_nt` / their specific BLAST subject DB.** We validate against kingdom-partitioned RefSeq organelle DBs (§4.2).
- **Their proposed-taxonomy panel.** Downstream of our boundary — Taxodactyl produces taxonomy; our report reports homology hits only.
- **Their fixed gene panel details.** Loci come from Taxodactyl's accepted-loci config at runtime ([brief.md §3.7](brief.md)).

### 6a.4 Run-level `run-report.html`

Uses the same visual language as the per-sample report but content is a
cross-sample overview:

- Sample count + status breakdown: `ok` / `low_coverage` (assembled below the warn floor — results are real but incomplete, [§2.1.3](02-stages.md#213-decision-matrix)) / `no_assembly` (recruited but nothing assembled) / `no_barcode` (assembled but no locus extractable) / `fail` (soft-failed at COVERAGE_GATE, no assembly attempted) / `error` (unexpected tool crash). This is the vocabulary of [CONSTITUTION principle 7](../CONSTITUTION.md), read from each sample's `metadata.json` `sample_status` field; `error` is the one value C7 derives itself, since a sample that crashed never produced a `metadata.json`. **`no_assembly` and `no_barcode` are separate states** and must not be merged into a single "no recovery" bucket: the first means nothing was assembled, the second means an organelle *was* assembled and is shipped as a real standalone output — only the barcode claim is negative. **`fail` and `error` are likewise not the same event**: `fail` is a coverage decision the pipeline made deliberately about the sample, `error` is the pipeline breaking. Label them so an operator can tell "this sample was too shallow" from "this run has a bug".
- Per-sample summary table — one row per sample_id, columns: kingdom, gate status, assembly outcome, coverage, top BLAST hit, panel recovery count. Each row links into `<sample_id>/report.html`. `low_coverage` rows are visually distinct from `ok` rows — a batch skim must not let a warned partial recovery read as a clean one.
- Run-level provenance panel: reference bundle version (§4.4 `manifest.json`), pipeline commit, samplesheet path + hash, invocation timestamp.

### 6a.5 Report implementation boilerplate

The renderer (task 43a_report_scaffold.md) is the boilerplate module,
extended rather than redesigned, split across three locations:

- **`bin/report/`** — the Python package. `report.py` exposes a single
  `render()` entry point, called in-process by `collate.py` (C6);
  `render_report.py` at `bin/`'s top level is the standalone CLI
  wrapper task 45_run_report.md (C7) reuses. Wrapped in a try/except at
  the call site — a rendering defect yields a minimal fallback report
  naming the sample, never a crash that takes the whole sample's
  bundle down (task 43a §5.2).
- **`scripts/report/templates/` + `scripts/report/static/`** — Jinja
  templates and the vendored Bootstrap 5 / Plotly / DataTables assets,
  staged into `COLLATE` (and `RUN_REPORT`) as explicit `path` inputs
  from a value channel in `main.nf`, not via `bin/`'s auto-staging —
  `bin/` is staged onto *every* process, and the 1.6 MB of vendored JS
  is only needed by the two report-emitting stages.
- **Self-contained single-file HTML** — all CSS, JS, and images inlined
  at render time. No external asset fetch; report is viewable offline
  and archivable as-is.
- **`metadata.json` is the only JSON input** (task 43a §3) — the prior
  glob-based `Config` result-file discovery, keyed off a per-sample
  `RESULT_DIR`, is retired. Task 42 already inlines every diagnostic
  payload into `metadata.json`; re-deriving values by globbing the work
  directory would give the report a second, unvalidated path to the
  same numbers (rule 18).
- **Class-based result objects** (`results.py`): `AbstractDataRow` for
  single-row typed structs and `AbstractResultRows` for tabular results
  driven by a CSV schema in `schema/`. Kept from the boilerplate for
  task 43b_report_stage_tabs.md's tabular stage content; `metadata.json`
  fields go straight into the render context, so no concrete
  `AbstractDataRow` subclass exists yet.
- **Tabbed layout (`content-tabs.html`)**, with a reusable
  analyst-comment macro (`subjective.html`) and a save-as-read-only
  workflow (`save-modal.html` + `save-report.js`). The strip carries
  §6a.2's four reader-perspective tabs — Overview (default), Assembly,
  Validation, Barcodes — **not** one tab per pipeline stage. Task 43a
  fixed the strip's shape against the earlier eight-stage layout; task
  43b_report_stage_tabs.md restructures it to the above and fills the
  panes.
- **Renderer inputs beyond `metadata.json`.** The file-shaped artefacts
  are passed as separate arguments and inlined at render time: the
  organelle map SVG (injected into the DOM, not via `<img>`, so
  per-feature tooltips bind), the Bandage graph PNG and annotation GFF
  (base64 `data:` URIs), the two NanoPlot HTML reports, and
  `barcodes.fasta` for §6a.2 tab 4's FASTA download. Wall-time
  reporting additionally needs the workflow start timestamp passed
  through `COLLATE`, since `metadata.json` carries no timing data.
- **`params.facility` and `params.analyst_name`** supply the submitting
  lab and analyst for the Overview tab's run-context panel. Both
  default to null and render as `-` when unset — an unset field must
  read as "not supplied", never as an empty row that looks like a
  missing result.


### 7 Specific feature requests

- Annotation GFF file should be available for download
