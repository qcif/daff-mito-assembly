# bin/

Helper scripts invoked by Nextflow processes. Nextflow auto-stages this
directory onto every process `PATH`, so scripts placed here are callable
by name inside any container.

One file per custom-logic component (plan.md §2.2):

| Script | Component | Stage | Phase |
|---|---|---|---|
| `parse_samplesheet.py` | C1 | 0 PARSE_SAMPLESHEET | P1 |
| `coverage_gate.py` | C2 | 6 COVERAGE_GATE | P1 |
| `bin_target.py` | C3 | 10 BIN_TARGET | P3 |
| `plastid_canonicalise.py` | C4 | 10 BIN_TARGET (plant-cp) | P3 |
| `validate_barcodes.py` | C5 | 13 MINIPROT_EXTRACT | P3 |
| `select_genetic_code.py` | C9 | 12 SELECT_GENETIC_CODE (animal_mt clade trial) | P4a |
| `collate.py` | C6 | 15 COLLATE | P4 |
| `report/` | C6/C7 | 15 COLLATE, 16 RUN_REPORT | P4a |
| `render_report.py` | C6/C7 | 15 COLLATE, 16 RUN_REPORT | P4a |
| `run_report.py` | C7 | 16 RUN_REPORT | P4 |

Requirements are documented in ../scripts/requirements.txt
