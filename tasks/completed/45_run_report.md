# Task 45 — Stage 16 `RUN_REPORT`: run manifest + cross-sample summary (C7)

**Phase:** P4a (from spec §6) — the last stage in
the pipeline, and the last of the P4a tail.
**Depends on:** task 42 for the
`metadata.json` it reads, task 43a_report_scaffold.md for the Jinja
machinery it reuses.
**Goal:** Replace the P0 stub in
`modules/local/run_report.nf` with real
`bin/run_report.py` (C7) emitting `run_manifest.json` and
`run-report.html`.

Smallest of the four. The join point across all samples — it does not
care whether each sample succeeded or soft-failed, only that a
`metadata.json` exists
(spec §2.1.3).

## Scope sketch

Per spec §6a.4:

- **Status classification** into `ok` / `low_coverage` / `no_assembly` /
  `no_barcode` / `fail` / `error`. Note this is **not** the vocabulary
  currently written in spec 06a-reports.md §6a.4 — that list says
  `no_recovery` and omits `no_assembly`. Task
  42_collate_bundle_metadata.md §3.3 establishes that `no_recovery` is
  stale 2026-07-27 text that survived a later rename sweep, and retires
  it; C7 reads the five-value `sample_status` from `metadata.json` and
  adds `error` for samples that crashed without producing one.
- **Per-sample summary table** — one row per sample_id: kingdom, gate
  status, assembly outcome, coverage, top BLAST hit, panel recovery
  count. Each row links into `<sample_id>/report.html`.
- **`run_manifest.json`** — samplesheet snapshot + hash, reference-bundle
  version (spec §4.4),
  pipeline commit, invocation timestamp.
- Same visual language as the per-sample report; reuse task 43a/43b's
  templates rather than forking them.
- Real container: shares `wf5/report` with C6 — replace the
  `python:3.12-slim` + `TODO P4` placeholder in
  `conf/containers.config`.

## Things this must not get wrong

- **`fail` and `error` are not the same event** and must never be merged
  in the display. `fail` is a coverage decision the pipeline made
  deliberately *about the sample*; `error` is the pipeline *breaking*.
  Label them so an operator can tell "this sample was too shallow" from
  "this run has a bug".
- **`low_coverage` rows are visually distinct from `ok` rows.** A batch
  skim must not let a warned partial recovery read as a clean one.

## Open questions for the expansion pass

- How is `error` actually detected? `COVERAGE_GATE` uses
  `errorStrategy 'ignore'` so failure is data, not a Nextflow error — but
  an unexpected crash in a *downstream* stage means no `metadata.json` is
  produced at all for that sample. Does C7 reconcile against the
  samplesheet to notice a sample that vanished, or does something upstream
  guarantee a bundle always exists?
- Acceptance for P4a includes rendering "a mixed batch with one
  soft-failed sample" — decide whether that becomes an integration
  fixture combination or a stub-profile test.

## Outcomes

**Open questions resolved:**

- **`error` detection: reconcile against the samplesheet.**
  `modules/local/collate.nf`'s own comment (task 42) already guarantees
  COLLATE emits a `metadata.json` for every sample it reaches, including
  the `fail`/`no_assembly` branches — so a samplesheet `sample_id` with
  no matching bundle in C7's collected set crashed somewhere *upstream*
  of COLLATE. `bin/run_report.py` parses the samplesheet CSV directly
  (it's staged into RUN_REPORT already) and builds one row per
  samplesheet `sample_id`, in samplesheet order; a `sample_id` absent
  from the collected `metadata.json` set — or present but unparseable —
  becomes `error`. A malformed/corrupt bundle is treated the same as a
  missing one: C7 cannot distinguish "crashed" from "produced garbage",
  and both mean no trustworthy record survived.
- **Mixed-batch acceptance: unit tests, not a new integration fixture.**
  `tests/integration/assertions.sh` (§ "Derived sample_status expected
  from COLLATE") already documents that no fixture reaches
  `fail`/`no_assembly`/`no_barcode` — a known, tracked gap from task 42
  §9, not something to backfill here. `scripts/tests/test_run_report.py`
  exercises all six `sample_status` values (five real + derived `error`)
  against synthetic `metadata.json` fixtures in one mixed batch
  (`TestBuildManifest.test_mixed_batch_reconciles_and_counts`), which is
  where P4a's "mixed batch with one soft-failed sample" acceptance is
  actually satisfied. The real integration fixtures (`INT-ANIMAL-01` ok,
  `INT-PLANT-01-pt` ok, `INT-PLANT-01-mt` low_coverage) already form a
  mixed `ok`/`low_coverage` batch and now carry their own
  `run_manifest.json`/`run-report.html` content assertions (added to
  `tests/integration/assertions.sh`: schema validation, a link to every
  sample's `report.html`, self-containment).

**Deviations / extensions beyond the stub's scope sketch:**

- **Reference-bundle manifest inlined in full.** Spec §4.4 states
  outright that "the full manifest is RUN_REPORT's (C7) to inline into
  `run_manifest.json`" (COLLATE copies only `version`/`generated_at`
  into each sample's `metadata.json`) — `RUN_REPORT` now takes
  `refs_manifest` as an explicit input (already a value channel in
  `main.nf`, reused rather than duplicated) and `run_manifest.json`
  embeds it verbatim under `reference_bundle`.
  `assets/run_manifest.schema.json` (new, mirroring
  `sample_metadata.schema.json`'s convention) is validated the same way
  COLLATE validates `metadata.json` — warn, never fail the run.
- **Resolved the `pipeline_commit: "unknown"` carry-forward item**
  (`tasks/todo.md`, filed by task 43a 2026-09-02, explicitly assigned to
  this task — now removed from `todo.md`). `main.nf` gains
  `resolvePipelineCommit()`: falls back from `workflow.commitId` /
  `workflow.revision` to a live `git rev-parse HEAD` against
  `workflow.projectDir` before giving up and returning `"unknown"`. This
  resolves a real commit hash for the everyday local-dev case (running
  from a working directory rather than a tagged checkout), which is what
  the per-sample provenance panel was silently eating before. Scoped to
  `RUN_REPORT` only, per the todo item's own framing ("task 45 owns
  *run-level* provenance") — `COLLATE`'s per-sample `metadata.json`
  still uses its original inline fallback
  (`workflow.commitId ?: workflow.revision ?: "unknown"}`) and was
  deliberately left untouched rather than opportunistically edited on a
  completed task's module.
- **`run_manifest.json` carries a full samplesheet snapshot, not just a
  hash.** `spec/00-overview.md`'s output tree comment and CONSTITUTION
  rule 16 both say "samplesheet snapshot"; a SHA256 alone doesn't let a
  reader six months from now see what was actually submitted if the CSV
  on disk has since moved or changed. `samplesheet.rows` (the parsed CSV,
  in samplesheet order) sits alongside `samplesheet.path` and
  `samplesheet.sha256`.
- **`conf/containers.config`'s `RUN_REPORT` entry needed no change.**
  The task brief's scope sketch says to "replace the `python:3.12-slim` +
  `TODO P4` placeholder" — that placeholder no longer exists; tasks 42/43a
  already pointed `RUN_REPORT` at the shared `neoformit/daff-wf5-scripts`
  image in advance (§6a.5: "C6/C7/the Jinja report renderer share one
  image"), anticipating this task. Verified via the CI lint job's
  container-coverage check, run locally.

**Design not spelled out in the stub, decided during implementation:**

- Reused `bin/report/`'s Jinja machinery (`get_static_file_contents` for
  self-contained-HTML asset inlining) via a new sibling module
  `bin/report/run_report.py`, rather than extending `report.py`'s
  `render()`/`build_context()` — the run-level page has a completely
  different content model (no tabs, no per-sample warnings mirror, no
  organism-specific title) and forcing it through the per-sample
  `build_context()` would have meant threading a lot of `None`-guarding
  through code that assumes a single sample's `metadata.json`. New
  templates: `scripts/report/templates/run_index.html` +
  `components/run-heading.html` / `run-summary.html` /
  `run-provenance.html`, styled with the same Bootstrap table classes
  (`table-success`/`table-warning`/`table-danger`) the per-sample report
  already uses for its Key findings table.
- `fail` and `error` share `table-danger` styling (both are genuinely bad
  outcomes) but always carry distinct label text ("Failed coverage gate"
  vs. "Pipeline error") — spec §6a.4's "must never be merged in the
  display" is read as a labelling requirement, not a colour requirement;
  Bootstrap's contextual palette doesn't have six distinguishable
  severities to spend on six statuses.
- `bin/run_report.py` mirrors `collate.py`'s exit-0/never-raise contract
  (CONSTITUTION rule 8): schema validation warns to stderr rather than
  failing, and a Jinja render failure falls back to a minimal HTML page
  naming the traceback while `run_manifest.json` — the audit-trail
  artefact — still gets written normally.

**Verification:**

- `bash scripts/pytest.sh` — 581 tests pass, 100% branch coverage
  project-wide (`bin/run_report.py` and `bin/report/run_report.py` both
  100%).
- `/home/cameron/.local/envs/claude/bin/flake8` clean on all new/changed
  Python.
- `nextflow run . -profile stub -stub-run` — the locally installed
  Nextflow 26.04.6 fails to *compile* `main.nf` on an unrelated
  pre-existing multi-line-string continuation at line 58
  (`validateParams()`), reproduced identically on a clean `git stash` of
  this branch's changes, so it is not something this task introduced —
  flagged for the user, not fixed here. Re-ran with `NXF_VER=25.10.2`
  (one of the other locally cached toolchains) and the full stub DAG —
  including `RUN_REPORT` with its new
  `refs_manifest`/`run_manifest_schema`/`pipeline_commit` inputs —
  completed green. Cleaned up with
  `.claude/scripts/clean_nextflow_run.sh` afterwards.
- Did not run `-profile integration` (real tools, real ONT fixtures,
  nightly cadence per rule 15) — the new `assertions.sh` block was
  syntax-checked (`bash -n`) and hand-verified against the same
  `run_manifest.json`/`run-report.html` shape the unit tests already
  exercise, but a real end-to-end run against the integration fixtures
  is worth doing before this ships to confirm the real-world render
  (left for the user to trigger via the `remote-run` skill or CI, given
  the shared-resource cost of an ad-hoc remote run).
