## 6a. Report design (per-sample `report.html`)

**Design reference:** [`mock-report.html`](mock-report.html). Use its **visual
layout, section structure, and interaction patterns** as the template.
**Do not** copy its scientific spec — the mock was produced for a SPAdes /
BUSCO-based pipeline, and its tool choices (assembler, completeness metric,
BLAST target DB, gene panel) must not leak into our specification. Content
comes from *our* pipeline stages (§2); the mock only informs how it's
presented.

### 6a.1 Design elements to adopt

- **Bootstrap-based single-file HTML** — self-contained, no external asset fetch, viewable offline.
- **Two-column header row:** *Key findings* (left, bulleted plain-English summary) + *Input parameters* (right, sample_id / kingdom / submitter metadata table).
- **Modal dialogs for detail views** — e.g. full parameter dump, tool versions with commit/build hashes. Keeps the main flow terse; detail one click away.
- **Multiple H2 top-level sections**, each a distinct part of the story, ordered as the pipeline runs.
- **Interactive charts** where a plot beats a table — read-length distribution, coverage-along-genome, per-contig coverage. Plotly is the mock's choice; keep unless it turns out to bloat the file substantially at multi-sample scale.
- **Inline SVG for the annotated organelle map** — clickable/hoverable gene features with tooltip metadata (name, strand, coordinates, source).
- **"Negative clarity" front-and-centre** ([brief.md §3.8](brief.md)) — when no organelle is assembled or no locus is extractable, the *Key findings* block says so explicitly at the top of the report, not buried in section 6.

### 6a.2 Section outline (our pipeline's content)

| H2 section | Content source | Notes |
|---|---|---|
| Header (H1 + two-column summary) | Sample metadata + top-line results | Key findings bullets: assembly outcome, coverage, top BLAST hit, gene-feature count, panel-locus extraction count. Right column: sample_id, kingdom, and any submitter-supplied optional columns from `samples.csv` (§0). |
| Sequencing data quality | `NANOPLOT_RAW` + `NANOPLOT_CLEAN` | Pre/post read stats, length + Q-score distributions, filter-yield summary from CHOPPER + FILTLONG. |
| Recruitment + coverage gate | `RECRUIT` + `COVERAGE_GATE` | Recruited-read count and total bases, estimated coverage vs configured MIN/MAX (§2.1), gate decision (pass / subsampled / soft-fail) with the seqtk fraction and seed if applicable. On plant targets, also the sibling-organelle split of the recruited pool from `sample_status.json` — `target_assigned_bases` vs `sibling_assigned_bases` and the resulting `sibling_organelle_fraction` — so the reader sees what the estimate was corrected for ([§2.1.5](02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)). Where `coverage_basis` is `total_recruited`, say so: the sibling split was unavailable and the estimate may over-state target depth. When soft-failed, this section is the terminal content — no assembly attempted. |
| *De novo* organelle assembly | `METAFLYE` | Assembly statistics (contigs, N50, total bp), per-contig breakdown (length, coverage, circularised y/n). Polishing is Flye's built-in single-iteration pass — no separate polish section until MEDAKA is un-deferred (§7 Q6). |
| Assembly quality assessment | `BANDAGE_NG` + `BIN_TARGET` | Assembly graph image, per-contig coverage-along-genome plot, binning classification (target / secondary / off-target). Carries the plastid QC signals of [§3.6](03-organelles.md#36-plastid-quadripartite-canonicalisation): the edge-count reading, the "canonical structure, no taxonomic support" warning when the `path1` substitution was withheld, and `target_source` provenance when it was applied. **No BUSCO** — mock includes it, our spec does not (kingdom-agnostic organelle BUSCO sets are patchy; revisit if needed). |
| Homology to reference databases | `BLAST_VALIDATE` | Top hits per target contig against the kingdom-appropriate RefSeq organelle DB (§4.2). Identity, coverage, subject accession. Explicit note when top-hit identity is below species-threshold — flags novel/underrepresented taxa without pretending to *assign* them (assignment is Taxodactyl's job — see [brief.md §2](brief.md) boundary). |
| Extracted barcode panel | `EXTRACT_BARCODES` | Per-locus: extracted length, ORF status, genetic code used, coordinates, protein-identity to the miniprot reference. Panel loci not recovered are listed with the reason (not present / ORF broken / below length threshold). |
| Annotated organelle genome map | `ANNOTATE` | Inline SVG. CDS features come from the broad miniprot pass (stage 12) for every target — the same features the barcode panel above is drawn from, so a locus shown here and a locus shipped in `barcodes.fasta` are the same object. tRNA/rRNA features come from MITOS2 on `animal_mt`; plant arms are CDS-only pending [§8 item 3](07-open-questions.md#8-remaining-open-questions) and **must say so** rather than leaving the absence to be inferred ([principle 7](../CONSTITUTION.md)). Tooltip shows name, product, strand, coordinates, feature source. Where `annotation_summary.json` records `cds_crosscheck` disagreements between miniprot and MITOS2, surface them — two independent methods disagreeing on a gene is a real QC signal. |

### 6a.3 Elements to leave in the mock

- **SPAdes references / SPAdes-specific stats.** We use metaFlye — swap all assembler mentions.
- **BUSCO section.** Not in our spec. Do not add.
- **`core_nt` / their specific BLAST subject DB.** We validate against kingdom-partitioned RefSeq organelle DBs (§4.2).
- **Their proposed-taxonomy panel.** Downstream of our boundary — Taxodactyl produces taxonomy; our report reports homology hits only.
- **Their fixed gene panel details.** Loci come from Taxodactyl's accepted-loci config at runtime ([brief.md §3.7](brief.md)).

### 6a.4 Run-level `run-report.html`

Uses the same visual language as the per-sample report but content is a
cross-sample overview:

- Sample count + status breakdown: `ok` / `no_recovery` (assembled but no locus) / `low_coverage` (soft-failed at COVERAGE_GATE) / `hard_failure` (unexpected tool crash).
- Per-sample summary table — one row per sample_id, columns: kingdom, gate status, assembly outcome, coverage, top BLAST hit, panel recovery count. Each row links into `<sample_id>/report.html`.
- Run-level provenance panel: reference bundle version (§4.4 `manifest.json`), pipeline commit, samplesheet path + hash, invocation timestamp.

### 6a.5 Report implementation boilerplate

The [`wf-report-boilerplate/`](../wf-report-boilerplate/) module is the
starting point for both report tiers. It is trimmed from a prior
Nextflow workflow to a clean scaffold that preserves the paradigms we
want to reuse:

- **Python + Jinja2 templating** with a single `render()` entry point
  ([`report.py`](../wf-report-boilerplate/report.py)) invoked by the
  `COLLATE` / `RUN_REPORT` processes.
- **Self-contained single-file HTML** — all CSS, JS, and images inlined
  at render time via `_get_static_file_contents()`. No external asset
  fetch; report is viewable offline and archivable as-is.
- **Config-driven result-file discovery** ([`config.py`](../wf-report-boilerplate/config.py)):
  glob-based `@property` accessors on `Config`, keyed off a per-sample
  `RESULT_DIR`. One property per output artefact; add as stages emit
  new files.
- **Class-based result objects** ([`results.py`](../wf-report-boilerplate/results.py)):
  `AbstractDataRow` for single-row typed structs (e.g. `Metadata`,
  `RunQC`) and `AbstractResultRows` for tabular results driven by a
  CSV schema in `schema/`. Keeps template logic dumb; casting,
  bootstrap-class assignment, and derived fields live on the class.
- **Bootstrap 5 + Plotly + DataTables** shipped as vendored static
  assets; tabbed layout (`content-tabs.html`), stage-per-tab, with a
  reusable analyst-comment macro (`subjective.html`) and a
  save-as-read-only workflow (`save-modal.html` + `save-report.js`).

The boilerplate ships stubbed — concrete result classes, config
accessors, and stage-per-tab component templates are added as each
pipeline stage (§2) lands. [`mock-report.html`](../wf-report-boilerplate/mock-report.html)
is retained inside the module as the visual design reference described
in §6a.
