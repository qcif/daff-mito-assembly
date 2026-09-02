# Task 43a — Stage 15 `COLLATE`: report renderer scaffold + Overview tab

**Phase:** P4a (spec 06-phases.md).
**Goal:** Turn `scripts/report/boilerplate/` — currently ~550 lines of
unused Python plus vendored Bootstrap 5, Plotly, DataTables and the DAFF
logo — into a working report renderer wired into `COLLATE`, and ship the
**header block + Overview tab**. Every input the full report will ever
need is plumbed here; the seven stage tabs are task
43b_report_stage_tabs.md.

**Depends on:** task 42_collate_bundle_metadata.md, which shipped the
`metadata.json` this renders from. That contract is fixed and
schema-validated (`assets/sample_metadata.schema.json`).

**Scope boundary — this task ships the *machinery and the frame*.**
It ends with a real, populated `report.html` carrying the always-visible
Key findings block, the Input parameters column, and the Overview tab
with its warning mirror. That is deliberately the part CONSTITUTION.md
principle 7 leans on hardest: negative clarity lives in Key findings, and
it should not wait behind seven tabs of stage detail.

**Related tasks:** 43b_report_stage_tabs.md builds the stage tabs on this
scaffold. 45_run_report.md reuses this task's Jinja machinery — it needs
the renderer, not the tabs, so it is unblocked by 43a alone.
44_organelle_map.md is independent; both report tasks build against its
current stub.

---

## 1. Overview

Every sample in this pipeline now ends with a directory containing an
assembled organelle, its annotation, the barcodes destined for
Taxodactyl, and a 20 kB `metadata.json` recording how all of it was
produced. What it does not have is a way for a person to read any of it.
`report.html` is currently a zero-byte file.

That file is the workflow's primary human-facing output. A biosecurity
officer opening it needs to answer, in order: did this sample work, how
confident should I be, and what exactly do I have. This task answers the
first question — properly, and before anything else — and builds the
machinery that lets 43b answer the other two.

The data layer is done. `metadata.json` already inlines the coverage
gate's whole decision, the binning metadata, per-contig assembly stats,
top BLAST hits, per-locus barcode outcomes, the full annotation summary,
and the provenance block. Task 42 built it first on purpose, so this task
is rendering, not re-derivation.

---

## 2. The boilerplate is the design reference — there is no mock

**Settle this first, because the spec is wrong about it.**
spec/06a-reports.md §6a opens by naming `mock-report.html` as the design
reference for the entire visual layer, and §6a.5 names it again. **That
file has never existed in this repository** —
`git log --all -- '*mock-report*'` returns nothing. Every §6a.5 path also
still points at `../wf-report-boilerplate/`, a location retired in commit
`50c60c3` ("Move report boilerplate") when the module moved to
`scripts/report/boilerplate/`.

**Decision (confirmed with the project owner, 2026-09-02):** the
boilerplate module *is* the design reference. It is trimmed from a prior
Nextflow workflow specifically to preserve the paradigms this project
should extend, and the point of extending it is that wf5's report is
structured like the other workflows in the family. There is no mock to
reconstruct and none is needed — §6a.1 and §6a.2 already enumerate the
layout in prose, and the boilerplate already implements most of the
chrome that prose describes.

Consequences:

- **Extend the boilerplate; do not redesign it.** Its tabbed layout, its
  `AbstractDataRow` / `AbstractResultRows` split, its
  `_get_static_file_contents()` inlining, its analyst-comment macro and
  its save-as-read-only workflow all stay. Fill in the stubs.
- §6a.3's "elements to leave in the mock" (SPAdes stats, BUSCO,
  `core_nt`, their proposed-taxonomy panel, their fixed gene panel) still
  stands as a **negative** constraint — those things must not appear —
  even though the mock they were written against is absent. Nothing in
  the boilerplate carries them, so this is a "do not add" list, not a
  "strip out" list.
- **Fix the spec** (§10).

---

## 3. Architecture: `metadata.json` is the only input

The boilerplate's `config.py` discovers result files by **globbing a
`RESULT_DIR`** — one `@property` per artefact, `sample_id` derived from
the name of a canonical result file. That model predates task 42. It is
now the wrong model, and following it would be the main way this task
goes wrong.

Task 42 §5.4 resolved to **inline** every diagnostic payload into
`metadata.json` rather than reference sibling files. So the report's
content is a single 20 kB JSON document that is already schema-validated,
already normalised, and already the Taxodactyl handoff record.
Re-deriving any of it by globbing the work directory would give the
report a second, unvalidated path to the same numbers — exactly the
divergence CONSTITUTION.md rule 18 exists to prevent.

**Decision:** the renderer takes `metadata.json` as its input contract.
`sample_id`, `kingdom`, `organelle`, `sample_status`, coverage,
assembly, homology, barcodes, annotation and provenance all come from
there and nowhere else.

Retire the glob layer: delete the `_get_file_by_pattern` result
accessors, `Config.sample_id`'s `NotImplementedError` stub,
`Config.DEFAULT_PARAMS_PATH`, `Config.VERSIONS_PATH`, `Config.FLAGS_CSV`
and `Config.flags` — all of which point at files that do not exist in
this repo and never will. An accessor that always returns `{}` is a
silent lie (rule 18).

---

## 4. Plumbing — every renderer input is wired in this task

**This is the task-42 lesson applied to itself:** "the renderer should
not be the thing that discovers its inputs are unreachable" (task 42
§2.2). 43a fixes every input defect and fixes the renderer's **complete
CLI surface**, including arguments 43b will be the first to actually
consume. The contract must not change shape twice.

### 4.1 NanoPlot is staged into `COLLATE` but never read

`main.nf` joins `NANOPLOT_RAW.out.reports` and
`NANOPLOT_CLEAN.out.reports` into `COLLATE`'s input tuple, and
`modules/local/collate.nf` declares both as `path(...)`. **`collate.py`
has no `--nanoplot-raw` or `--nanoplot-clean` argument.** The
directories are staged into the task sandbox on every sample and then
ignored, so spec §6a.2's "Sequencing data quality" section has no data
source at all.

`NanoStats.txt` (NanoPlot runs with `--tsv_stats`) is a two-column TSV
and parses trivially. A real one:

```
number_of_reads	1242
number_of_bases	7900845.0
median_read_length	3374.0
mean_read_length	6361.4
read_length_stdev	7065.2
n50	13081.0
mean_qual	11.0
median_qual	11.5
Reads >Q10:	1002 (80.7%) 6.1Mb
```

Add `--nanoplot-raw` / `--nanoplot-clean` to `collate.py`, parse both
`NanoStats.txt` files, and add a **`read_qc`** block to `metadata.json`
carrying `raw` and `clean` sub-objects plus a derived `filter_yield`
(reads and bases retained, as counts and percentages). Extend
`assets/sample_metadata.schema.json` accordingly — nullable, since a
sample that soft-fails at the gate still has both NanoPlot runs but a
malformed one must not break collation.

Deriving filter yield from raw-vs-clean read and base counts is what
satisfies §6a.2's "filter-yield summary from CHOPPER + FILTLONG" without
adding a stats-emitting step to either tool. Record that reasoning.

**Do not re-render NanoPlot's distribution plots.** NanoPlot emits them
as standalone Plotly HTML files (`WeightedHistogramReadlength.html` and
friends), each embedding its own copy of plotly.js. Scraping trace JSON
out of them would be fragile against every NanoPlot release — rule 19's
drift problem. Instead, pass `NanoPlot-report.html` (~192 kB) to the
renderer for embedding as a base64 `data:` URI; 43b puts it behind a
"full read QC report" button.

### 4.2 `recruit_stats.json` never reaches `COLLATE`

`RECRUIT` publishes `recruit_stats.json` (`reads_aligned`,
`reads_recruited`, `min_aligned_frac`, `min_aligned_bp`) but it is not in
`COLLATE`'s input tuple. §6a.2's recruitment section asks for
"recruited-read count and total bases"; `metadata.json` currently has
only the bases, via `coverage.gate.total_recruited_bases`.

Wire it — one `.join()` in `main.nf`, one `path` in `collate.nf`, one
argument in `collate.py` — and put the recruitment thresholds actually
applied into the audit trail (rule 16). If it complicates the
`ch_failed_inputs` placeholder list disproportionately, dropping the read
count and reporting bases only is acceptable, but **record the choice**
rather than leaving 43b's section silently short of spec.

### 4.3 Renderer arguments for the file-shaped artefacts

Three artefacts are files, not JSON fields. Wire all three as optional
renderer arguments now; 43b renders them.

- **Organelle map SVG** — must be inlined into the DOM (not an `<img>`)
  so task 44's per-feature tooltips work. Today it is a zero-byte stub;
  the renderer must treat that as "not yet available" without erroring.
- **Bandage graph PNG** — base64 `data:` URI.
- **Annotation GFF** — base64 `data:` URI on a `download`-attributed
  anchor. This satisfies spec §6a "Specific feature requests" *and* keeps
  the report self-contained: a relative link to
  `../organelle_annotation.gff` breaks the moment the report is emailed
  on its own, which is exactly how these get shared.

### 4.4 Parameters modal

`inputs.html` already renders a "View all parameters" modal from a
`parameters` context dict and a "View tool versions" modal from
`versions`. Tool versions come free from
`metadata.json.provenance.tool_versions`. Parameters have no source.

Emit the resolved params from the module using the heredoc pattern
`collate.nf` already uses for `meta.json` — the same command-injection
reasoning applies, since params carry arbitrary paths and strings.

---

## 5. Packaging and failure containment

### 5.1 Where the renderer lives

The boilerplate is a **Python package** (`from . import config`), not a
flat script, and CONSTITUTION.md rule 12 forbids baking source into the
image. Two facts settle it:

1. **Nextflow stages `bin/` subdirectories, and a package under `bin/`
   imports cleanly.** Verified during this brief with a standalone
   Nextflow 25.10.2 repro: `bin/pkg/__init__.py` imported successfully
   from `bin/probe.py` via
   `sys.path.insert(0, str(Path(__file__).resolve().parent))`.
2. **The static assets should not live in `bin/`.** `static/js/` alone is
   1.6 MB of vendored jQuery, Bootstrap, Plotly and DataTables. Nextflow
   stages the *entire* `bin/` directory for *every* process on a remote
   executor, so putting it there would upload 1.6 MB per task across the
   whole run on Azure Batch — for the benefit of two processes.

**Decision:**

- Move the Python package to `bin/report/` and add a thin
  `bin/render_report.py` CLI entry point. **`chmod +x` it** — task 42's
  outcomes note records that a missing executable bit surfaced as exit
  126 from `COLLATE`, because Nextflow preserves source permissions
  rather than forcing them.
- Keep `templates/` and `static/` outside `bin/`, staged into `COLLATE`
  and `RUN_REPORT` as an explicit `path` input from a value channel in
  `main.nf`. That is rule 13's pattern verbatim, it works on every
  executor, and it charges the 1.6 MB to the two processes that need it.
  `report.py`'s `TEMPLATE_DIR` and `STATIC_DIR` become CLI arguments
  instead of module-level constants.
- Update `bin/README.md`'s script table, and CI's flake8 exclude path in
  `.github/workflows/tests.yml` (currently
  `--exclude=scripts/report/boilerplate/static/`). **Verified:** the
  boilerplate is already flake8-clean at the default 79 columns, so no
  re-wrapping is needed.

### 5.2 The renderer must never break the bundle

`collate.py` is contractually always-exit-0 (task 42 §7) because
CONSTITUTION.md principle 8 forbids one sample's failure from taking its
siblings down. Adding a Jinja render introduces a large new class of
exceptions — a template typo, a missing key, a malformed SVG.

**Wrap the render in `try/except`.** On failure, write a minimal fallback
`report.html` that states plainly that report rendering failed, names the
sample, and includes the traceback, then continue and exit 0. The bundle,
`metadata.json` and the Taxodactyl handoff all survive; the failure is
loud and visible in the artefact where an operator will actually see it
(rule 18), rather than swallowed or escalated into a run abort.

---

## 6. Layout: keep the tabs, mirror every warning into Overview

**Decision (confirmed with the project owner, 2026-09-02).** spec §6a.1
asks for "multiple H2 top-level sections", but the boilerplate ships a
tabbed layout (`content-tabs.html`, stage-per-tab) and §2's decision —
extend the boilerplate, stay structurally consistent with the sibling
workflows — governs. **Keep the tabs.**

The risk tabs carry is the one principle 7 cares about: a warning behind
a tab the reader never opens is a warning that did not happen. Two
mechanisms close it, and **both are built in this task** because 43b's
tabs depend on them existing:

1. **The header block is outside the tab strip and always visible.**
   `heading.html` renders above `content-tabs.html`. `overview.html`'s
   existing two-column row (`inputs.html` | `qc.html`) becomes §6a.1's
   two-column header row: **Key findings** on the left — the
   plain-English bulleted summary, including every negative state — and
   **Input parameters** on the right, which `inputs.html` already
   implements against `sample_id` / kingdom / the submitter-supplied
   optional columns. `qc.html`'s placeholder content is removed; its tab
   is 43b's.
2. **Overview is the first tab and the canonical warning mirror.**
   Establish the mechanism as a **`warnings` list in the render
   context** — each entry carrying its text, severity and owning tab —
   which Overview iterates. 43b appends to that list rather than
   inventing per-tab warning markup. A stage tab keeps its warning in
   context, but no warning exists *only* behind a tab.

Register the tab strip with all eight tabs in pipeline order (per §6a.2),
with tabs 2–8 as empty placeholders naming the stage they will carry.
That fixes the navigation contract once, so 43b fills panes rather than
restructuring the strip.

---

## 7. Negative clarity is in Key findings, at the top

This is what the task is really for.

`no_assembly` and `no_barcode` are stated in the Key findings block, in
plain English, above the tab strip — not buried in the tab that owns the
stage where the negative arose (brief.md §3.8, principle 7). `fail`
likewise: the sample never reached assembly, and Key findings says so.

**`no_assembly` and `no_barcode` must not read as the same outcome.**
`no_barcode` means an organelle *was* assembled and is shipped as a real
standalone output (principle 6) — only the barcode claim is negative.
`no_assembly` means nothing was assembled. A report that renders them
with the same words has merged two of principle 7's three signals.

**`low_coverage` is a warned partial result, not a failure.** Key
findings states it as a real recovery carrying a caveat, never as a
degraded or failed one (principle 7's second half). The detailed
warn-floor warning is 43b's tab 3, but the top-line framing is set here,
and it must not read as a negative.

Key findings also carries the top-line summary §6a.2 specifies: assembly
outcome, coverage, top BLAST hit, gene-feature count, panel-locus
extraction count — all available from `metadata.json` today.

---

## 8. Container

`conf/containers.config` already pins `COLLATE` **and** `RUN_REPORT` to
`neoformit/daff-wf5-scripts:684d981@sha256:2fcdfbdeac14630f0dce95f0941043e67e8229586559cefec911fee3449fd3b2`
— the shared image carrying Python + `jinja2` + `jsonschema`, deps only,
source staged at runtime (rule 12). **Verified: no image change is
needed.** If a new dependency turns out to be genuinely necessary,
justify it against rule 19 first — the boilerplate deliberately vendors
its front-end assets rather than pulling a chart library at runtime, and
the same restraint applies to the Python side.

---

## 9. Tests — `scripts/tests/test_report.py`

The renderer is custom logic under `bin/`, so rule 14 applies: **100%
branch coverage**, run via `scripts/pytest.sh`. Mock at the tool boundary
only (task 27) — build real `metadata.json` fixtures, render real HTML,
assert on the real output string.

- **All five `sample_status` values render** without raising, and each
  produces a Key findings block naming its own state.
- **`no_assembly` and `no_barcode` produce different text** (§7) — the
  principle-7 regression that matters most in this task.
- **`low_coverage` Key findings does not use failure wording** and is
  distinguishable from `ok`.
- **Self-containment:** the rendered HTML contains no `src="http`,
  `href="http` or `url(http` outside deliberate informational links.
  This is the assertion that stops a future edit from quietly
  reintroducing a CDN dependency.
- **`read_qc` parsing** — well-formed, malformed and absent
  `NanoStats.txt` all handled; `filter_yield` maths correct including
  the zero-reads edge case.
- **Missing optional artefacts** — zero-byte SVG, absent PNG, absent
  GFF, absent `read_qc` — render without raising.
- **Render failure produces the fallback report and still exits 0**
  (§5.2), asserted through `collate.py`.
- **Emitted `metadata.json` still validates** against the extended
  schema in every branch; `test_collate.py` stays green.

---

## 10. Integration reconciliation

- **Stub run.** `-stub-run` must stay green with the new
  templates/static `path` input and the widened `COLLATE` tuple. Add the
  new channel to `conf/stub.config` and `conf/integration.config`.
- **`tests/integration/assertions.sh`.** Task 42 added a per-sample
  block; extend it with: `report.html` is non-empty, contains the
  `sample_id`, contains the sample's `sample_status`, and contains no
  external asset references.
- **Assert on structure and status, not on numbers** that shift with
  reference-data rebuilds (rule 19).
- **Known gap, carried forward from task 42:** no fixture exercises
  `fail`, `no_assembly` or `no_barcode` with real data. §9's unit tests
  are the only coverage of those three Key-findings paths. Record this;
  do not paper over it.
- **`pipeline_commit` is `"unknown"`** in current real output —
  `workflow.commitId` is null when the pipeline runs from a working
  directory rather than a cloned project. The provenance panel must
  render that honestly rather than hiding the row; add a `tasks/todo.md`
  entry for task 45, which owns run-level provenance.

---

## 11. Spec updates

- **spec/06a-reports.md §6a** — remove the `mock-report.html` design
  reference and §6a.5's mock sentence (§2); the boilerplate is the
  reference. Keep §6a.3 as a "do not add" list and say so.
- **spec/06a-reports.md §6a.5** — repoint every
  `../wf-report-boilerplate/` path at `bin/report/` and the staged asset
  dir per §5.1.
- **spec/06a-reports.md §6a.1** — record that the tabbed layout is the
  chosen realisation of "multiple top-level sections", with §6's
  always-visible-header + Overview-mirror rule stated as the principle-7
  guarantee.
- **spec/02-stages.md §2.2** — repoint the C6 row's report-renderer link;
  note the `read_qc` addition to `metadata.json` in stage 15's row.
- **tasks/todo.md** — add the `pipeline_commit` entry from §10. Leave the
  task 24, task 31, `ORGANELLE_MAP` and `RUN_REPORT` entries alone; the
  first two are resolved by 43b, not here.

Per the task-authoring convention: when editing spec files, reference
tasks as plain text (`task 43a_report_scaffold.md`), never as markdown
links, and never with line anchors.

---

## 12. Acceptance criteria

1. `COLLATE` writes a populated `report.html` for all three integration
   fixtures; `bin/render_report.py` and `bin/report/` are flake8-clean at
   79 columns with 100% branch coverage under `scripts/pytest.sh`.
2. The report is a **single self-contained file** — no external fetch,
   asserted by test (§9), viewable offline.
3. The always-visible Key findings header renders every `sample_status`
   distinctly, with `no_assembly` ≠ `no_barcode` and `low_coverage` not
   framed as a failure (§7).
4. Input parameters column, parameters modal and tool-versions modal all
   render from `metadata.json` + the emitted params JSON (§4.4).
5. The eight-tab strip is registered in pipeline order; Overview renders
   the `warnings` mirror mechanism 43b will append to (§6).
6. `read_qc` is in `metadata.json` and in the schema; NanoPlot is no
   longer staged-and-ignored (§4.1).
7. The renderer's **complete** CLI surface is fixed, including the
   arguments only 43b consumes (§4.3) — the contract does not change
   shape twice.
8. A render failure yields the fallback report and `COLLATE` still exits
   0 (§5.2).
9. `-stub-run` green; `assertions.sh` extended and passing on the nightly
   integration profile.
10. No dangling `mock-report.html` or `wf-report-boilerplate/` references
    remain anywhere in the repo.

---

## 13. Out of scope

- **The seven stage tabs and all their content rules** — task
  43b_report_stage_tabs.md. This task registers the tabs and leaves the
  panes as named placeholders.
- **Charts and the Plotly size decision** — 43b. The vendored assets are
  inlined here; whether to keep Plotly is measured there.
- **The real organelle map SVG** — task 44_organelle_map.md. Build
  against the current zero-byte stub.
- **`run-report.html` / `run_manifest.json`** — task 45_run_report.md.
  Factor shared components (heading, modals, status vocabulary,
  contextual-class helpers) into `templates/components/` and `results.py`
  so task 45 reuses rather than forks them, but do not build the run tier
  here.
- **Sourcing fixtures for `fail` / `no_assembly` / `no_barcode`** —
  tasks/todo.md, unchanged; the gap is recorded in §10, not closed here.

---

## 14. Outcomes

Implemented as specified, with a few deviations recorded below.

- **Package layout.** `bin/report/` (package: `config.py`, `report.py`,
  `results.py`, `utils.py`, `filters/`) + `bin/render_report.py` (thin
  CLI wrapper) + `scripts/report/templates/` + `scripts/report/static/`
  + `scripts/report/schema/`, per §5.1. `chmod +x` applied to every
  `bin/*.py` file including the new package, matching repo convention
  (not strictly required for files that are only imported, but
  consistent with every existing `bin/*.py`).
- **`config.py`** trimmed to two string constants (`REPORT_TITLE`,
  `REPORT_SUBTITLE_HTML`) rather than kept as a class — there was no
  state left to hang a `Config` instance off once the glob-based
  accessors were retired (§3), so the class wrapper would have been
  pure ceremony.
- **`results.py`** keeps `AbstractDataRow`/`AbstractResultRows`/`FLAGS`
  for 43b; the `Metadata`/`RunQC` concrete stubs are deleted rather
  than reshaped, since `metadata.json` fields go straight into the
  render context (no per-row casting needed for scalars already typed
  by the schema).
- **Analyst/facility/wall-time tracking** (the old boilerplate's
  `walltime.html`, tied to a `*_start_timestamp.txt` file that has
  never existed in this repo) is dropped rather than ported — nothing
  in spec §6a or this brief asks for it, and carrying forward a dead
  file-glob would repeat the exact defect §3 retires elsewhere.
  `templates/components/walltime.html` and `qc.html` (superseded by
  the promoted two-column header row, §6) are deleted. Overview's
  header column is now `provenance.html` (pipeline commit, reference
  bundle version/build date, bundle kind) instead — real content
  available today, where wall-time tracking was not.
- **Parameters modal** shows a flat `{name: value}` table rather than
  the old default/user split: the workflow has one resolved parameter
  set per run, not a separate defaults-file and user-overrides-file to
  diff. `collate.nf` emits the whole resolved `params` map via the
  same heredoc-to-file pattern as `meta.json` (command-injection
  reasoning applies identically — params values are not
  shell-interpolated).
- **`report_templates`/`report_static` are plain `file("${projectDir}/...")`
  value channels in `main.nf`, not new `params.*` entries.** They are
  part of the pipeline's own code, not swappable reference data, so a
  params knob would be a knob with no reason to turn — no
  `conf/stub.config`/`conf/integration.config` changes were needed for
  them as a result (unlike `refs_manifest` etc., they resolve to the
  same repo-relative path regardless of profile). This is the one
  deviation from §5.1/§10's suggestion to add config entries for the
  new channels.
- **`RECRUIT.out.stats`** was missing its `meta` tag entirely (a
  pre-existing gap — the emit was `path("recruit_stats.json")` with no
  `val(meta)`, making it unjoinable and, correspondingly, never
  actually joined into anything). Fixed to `tuple val(meta),
  path("recruit_stats.json")` so §4.2's join works; both the read
  count and bases are carried through (no need to fall back to
  bases-only).
- **Warnings mirror seeded with real content**, not left as a bare
  mechanism: low_coverage, `coverage_basis == total_recruited`,
  `annotator_failed`, `genetic_code_agreement == false`, and
  `cds_crosscheck.annotator_only` all populate it today. The last two
  close the tasks/todo.md task-31 carry-forward item (surfacing
  `cds_crosscheck` disagreements / genetic-code disagreement) at the
  point the mechanism was built, rather than deferring it to 43b.
- **Self-containment test allowlist** needed a second entry beyond the
  repo subtitle link: the vendored `plotly-basic` JS ships a literal
  `href="https://plotly.com/"` string (Plotly's own modebar
  attribution), inert until a chart actually calls `Plotly.newPlot`
  with `showlink` — pre-existing in the vendored asset, not introduced
  by this task. Both `scripts/tests/test_report.py` and
  `tests/integration/assertions.sh`'s self-containment check allowlist
  it explicitly.
- **Verification:** `bash scripts/pytest.sh scripts/tests/` — 398
  passed, 100% branch coverage across `bin/*.py` and `bin/report/**`
  (coverage script's `--include` glob widened to catch the new
  subpackage). `flake8` clean at 79 columns. `nextflow run main.nf
  -profile stub -stub-run` green end-to-end including the widened
  `COLLATE` tuple. A real `-profile integration` run of `INT-ANIMAL-01`
  (single-sample samplesheet, real Docker containers, real reference
  bundle) completed in 11m42s and produced a genuine (non-fallback)
  1.83 MB `report.html` with `sample_status: ok`, populated `read_qc`
  (72.79% reads / 72.92% bases retained through CHOPPER+FILTLONG) and
  `coverage.recruitment`, verified by hand for self-containment,
  `sample_id`, and `sample_status` presence — the same checks now
  encoded in `tests/integration/assertions.sh`'s new report.html block
  for the nightly profile.
- **Known gap, unchanged from task 42:** no fixture exercises `fail`,
  `no_assembly` or `no_barcode` with real data (§10) — still only unit
  tests reach those three Key-findings paths against realistic data
  shapes.
