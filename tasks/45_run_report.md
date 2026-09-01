# Task 45 — Stage 16 `RUN_REPORT`: run manifest + cross-sample summary (C7)

> **STUB — expand before executing.**

**Phase:** P4a (from spec §6) — the last stage in
the pipeline, and the last of the P4a tail.
**Depends on:** task 42 for the
`metadata.json` it reads, task 43 for the Jinja
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

- **Status classification** into `ok` / `low_coverage` / `no_recovery` /
  `fail` / `error`.
- **Per-sample summary table** — one row per sample_id: kingdom, gate
  status, assembly outcome, coverage, top BLAST hit, panel recovery
  count. Each row links into `<sample_id>/report.html`.
- **`run_manifest.json`** — samplesheet snapshot + hash, reference-bundle
  version (spec §4.4),
  pipeline commit, invocation timestamp.
- Same visual language as the per-sample report; reuse task 43's
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
