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
- **Two-column header row:** *Key findings* (left, bulleted plain-English summary) + *Input parameters* (right, sample_id / kingdom / submitter metadata table). Always visible, above the tab strip — not inside a tab (task 43a §6).
- **Modal dialogs for detail views** — e.g. full parameter dump, tool versions with commit/build hashes. Keeps the main flow terse; detail one click away.
- **A tabbed layout, stage-per-tab**, ordered as the pipeline runs — the chosen realisation of "multiple top-level sections" (task 43a §6), kept from the boilerplate rather than flattened into H2 sections. Overview is tab 1 and mirrors every warning raised on any other tab, via a `warnings` list in the render context — so a caveat is never visible only behind a tab a reader might not open.
- **Interactive charts** where a plot beats a table — read-length distribution, coverage-along-genome, per-contig coverage. Plotly is the mock's choice; keep unless it turns out to bloat the file substantially at multi-sample scale.
- **Inline SVG for the annotated organelle map** — clickable/hoverable gene features with tooltip metadata (name, strand, coordinates, source, confidence score).
- **"Negative clarity" front-and-centre** ([brief.md §3.8](brief.md)) — when no organelle is assembled or no locus is extractable, the *Key findings* block says so explicitly at the top of the report, not buried in section 6. `no_assembly` and `no_barcode` must never read as the same outcome, and `low_coverage` must never read as a failure (task 43a §7).

### 6a.2 Section outline (our pipeline's content)

| H2 section | Content source | Notes |
|---|---|---|
| Header (H1 + two-column summary) | Sample metadata + top-line results | Key findings bullets: assembly outcome, coverage, top BLAST hit, gene-feature count, panel-locus extraction count. Right column: sample_id, kingdom, and any submitter-supplied optional columns from `samples.csv` (§0). |
| Sequencing data quality | `NANOPLOT_RAW` + `NANOPLOT_CLEAN` | Pre/post read stats, length + Q-score distributions, filter-yield summary from CHOPPER + FILTLONG. |
| Recruitment + coverage gate | `RECRUIT` + `COVERAGE_GATE` | Recruited-read count and total bases, estimated coverage vs the configured HARD MIN / WARN / MAX (§2.1), gate decision (pass / pass-with-warning / subsampled / soft-fail) with the seqtk fraction and seed if applicable. On a **pass with warning** (`low_coverage`, [§2.1.3](02-stages.md#213-decision-matrix)) this section carries a prominent warning that everything below it was produced under the warn floor and that missing genes and fragmented contigs are expected — the reader must not read incompleteness as a negative finding. State the estimate's optimism explicitly ([§2.1.6](02-stages.md#216-provenance-of-the-two-floors)): actual depth on the assembled contig is typically well below `estimated_cov`. On plant targets, also the sibling-organelle split of the recruited pool from `sample_status.json` — `target_assigned_bases` vs `sibling_assigned_bases` and the resulting `sibling_organelle_fraction` — so the reader sees what the estimate was corrected for ([§2.1.5](02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)). Where `coverage_basis` is `total_recruited`, say so: the sibling split was unavailable and the estimate may over-state target depth. When soft-failed, this section is the terminal content — no assembly attempted. |
| *De novo* organelle assembly | `METAFLYE` | Assembly statistics (contigs, N50, total bp), per-contig breakdown (length, coverage, circularised y/n). Polishing is Flye's built-in single-iteration pass — no separate polish section until MEDAKA is un-deferred (§7 Q6). |
| Assembly quality assessment | `BANDAGE_NG` + `BIN_TARGET` | Assembly graph image, per-contig coverage-along-genome plot, binning classification (target / secondary / off-target). Carries the plastid QC signals of [§3.6](03-organelles.md#36-plastid-quadripartite-canonicalisation): the edge-count reading, the "canonical structure, no taxonomic support" warning when the `path1` substitution was withheld, and `target_source` provenance when it was applied. **No BUSCO** — mock includes it, our spec does not (kingdom-agnostic organelle BUSCO sets are patchy; revisit if needed). |
| Homology to reference databases | `BLAST_VALIDATE` | Top hits per target contig against the kingdom-appropriate RefSeq organelle DB (§4.2). Identity, coverage, subject accession. Explicit note when top-hit identity is below species-threshold — flags novel/underrepresented taxa without pretending to *assign* them (assignment is Taxodactyl's job — see [brief.md §2](brief.md) boundary). |
| Extracted barcode panel | `EXTRACT_BARCODES` | Per-locus: extracted length, ORF status, genetic code used, coordinates, protein-identity to the miniprot reference. Panel loci not recovered are listed with the reason (not present / ORF broken / below length threshold). |
| Annotated organelle genome map | `ANNOTATE` → `JOIN_ANNOTATION_SCORES` | Inline SVG. CDS features come from the broad miniprot pass (stage 12) for every target — the same features the barcode panel above is drawn from, so a locus shown here and a locus shipped in `barcodes.fasta` are the same object. tRNA/rRNA features come from MITOS2 on `animal_mt`; plant arms are CDS-only pending [§8 item 3](07-open-questions.md#8-remaining-open-questions) and **must say so** rather than leaving the absence to be inferred ([principle 7](../CONSTITUTION.md)). Tooltip shows name, product, strand, coordinates, feature source, and confidence: `pident`/`qcovhsp`/`bitscore` from `annotation_summary.json`'s `cds_scores` (task 40), read as two independent axes — high identity + low coverage is a truncated fragment of a well-referenced gene, low identity + high coverage is a complete gene in an under-referenced clade, and the second must not read as failure. A feature with a `null` triplet (no blastp hit against its own gene's panel) renders as explicitly unscored, not as zero confidence. The report may sort or style features by score — that is the only consequence a score is permitted to have (task 40 §5); it never changes which features appear or a sample's status. Where `annotation_summary.json` records `cds_crosscheck` disagreements between miniprot and MITOS2, surface them — two independent methods disagreeing on a gene is a real QC signal. `status: "annotator_failed"` (task 41) must render distinctly from `ok_cds_only`: both ship a CDS-only annotation, but `ok_cds_only` is the honest provisional state of a target with no annotator configured, while `annotator_failed` is a configured annotator that ran and returned nothing — the `reason` field names the tool and exit code and must be shown, not just the CDS-only feature set. A reader must never mistake a tool failure for a provisional-by-design state. |

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
- **Tabbed layout (`content-tabs.html`), stage-per-tab**, with a
  reusable analyst-comment macro (`subjective.html`) and a
  save-as-read-only workflow (`save-modal.html` + `save-report.js`).
  Registered in pipeline order — Overview + the seven §6a.2 stage
  sections — task 43a fixes the strip's shape once; task
  43b_report_stage_tabs.md fills the panes.


### 7 Specific feature requests

- Annotation GFF file should be available for download
