## 6. Development phases

| Phase | Goal | Exit criteria |
|-------|------|---------------|
| **P0 — Scaffold** | Nextflow DSL2 skeleton, container plan, params schema, CI lint. | `nextflow run main.nf -profile test` completes a no-op end-to-end. |
| **P1 — Read prep + recruit + coverage gate** | Stages 1–6 wired with real test data. | Recruitment yield verified on plant test data at skim depth; coverage gate correctly subsamples an over-covered sample and soft-fails an under-covered one without blocking sibling samples. |
| **P2 — Assemble + polish + viz** | Stages 7–9. | Plant chloroplast assembled and visualised end-to-end. |
| **P3 — Bin + validate + extract** | Stages 10, 11, 13. | Barcode FASTA produced for plant + animal clean datasets; ORF validation passes. |
| **P4 — Annotate + collate + report** | Stages 12, 14, 15, 16; full output bundle, per-sample `report.html`, run-level `run-report.html`. | Taxodactyl handoff bundle validates against schema; full annotation produced; both report tiers render on multi-sample test invocation, including a mixed batch with one soft-failed sample. |
| **P5 — Contamination tolerance** | Minor-contamination datasets. | Acceptance criteria in §5 met. |
| **P6 — Negative-result clarity** | No-assembly and no-locus paths explicit in outputs. | Degraded-sample test yields explicit no-recovery diagnostic, distinct from null output ([brief.md §3.8](brief.md)). |
