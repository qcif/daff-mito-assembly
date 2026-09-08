# Organelle assembly workflow

This Nextflow workflow is being developed for DAFF Biosecurity for use taxonomic identification. 
The workflow assembles mitochondrial or plastid genomes from Nanopore skim-sequencing reads.

## Debugging the report renderer

`bin/render_report.py` can be run/debugged directly in VS Code against a real
sample's output, via the `render_report.py (frozen fixture)` config in
`.vscode/launch.json`. It points at the `INT-PLANT-01-pt` sample under
`tests/integration/output/`, so that fixture must exist first — run the
Tier 2 integration test to populate it (see `docs/testing.md`):

```bash
bash tests/integration/fetch_fixtures.sh
bash scripts/fetch_refs.sh v2026.08
nextflow run . -profile integration
```

The launch config runs with `.vscode/venv`, a lightweight local venv
(gitignored) containing just `jinja2` — `render_report.py` doesn't need the
full `scripts/requirements.txt` stack. Create it once with:

```bash
python3 -m venv .vscode/venv
.vscode/venv/bin/pip install jinja2==3.1.6
```

Rendered output goes to `tests/manual/output/render_report_debug/report.html`
(also gitignored).
