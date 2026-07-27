# Development Plan: Organelle Genome Assembly + Barcode Recovery

**Status:** Draft v0.2 · **Companion:** [brief.md](brief.md) · **Format:** Nextflow DSL2

The plan lives under [`spec/`](spec/), split by top-level section. This
page is a lightweight table of contents; anchors used by other repo
files (`plan.md#…`) resolve to the section headers listed below.

## Contents

- [spec/00-overview.md](spec/00-overview.md) — Overview + input sample sheet
- [spec/01-pipeline-flow.md](spec/01-pipeline-flow.md) — Pipeline flow + engineering constraints
- [spec/02-stages.md](spec/02-stages.md) — Stage detail, coverage gate, custom logic
- [spec/03-organelles.md](spec/03-organelles.md) — Organelle considerations (mt vs cp)
- [spec/04-reference-data.md](spec/04-reference-data.md) — Reference data
- [spec/05-test-data.md](spec/05-test-data.md) — Test data, tests, CI
- [spec/06-phases.md](spec/06-phases.md) — Development phases
- [spec/06a-reports.md](spec/06a-reports.md) — Report design
- [spec/07-open-questions.md](spec/07-open-questions.md) — Reconciliation + open questions + fine tuning

## Related

- [brief.md](brief.md) — original scoping brief
- [tasks/](tasks/) — per-task specs
- [spec/brief.md](spec/brief.md) — brief copy alongside the plan
