"""Render the cross-sample `run-report.html` from `run_manifest.json`
(C7, spec §6a.4, task 45).

Reuses the Jinja machinery `report.py` (task 43a) established for the
self-contained-HTML asset inlining (`get_static_file_contents`) rather
than forking it, but renders a different template set: the per-sample
report's four reader-perspective tabs answer "did this one sample
work"; this page answers "how did the batch do" — no tabs, just a
status breakdown, a per-sample table linking into each `<sample_id>/
report.html`, and a run-level provenance panel.
"""

from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from .config import REPORT_SUBTITLE_HTML
from .report import get_static_file_contents

# Five-value vocabulary of CONSTITUTION.md principle 7, plus `error`
# (spec §6a.4) — the value run_report.py derives itself for a
# samplesheet sample_id with no matching metadata.json.
STATUS_LABELS = {
    "ok": "OK",
    "low_coverage": "Low coverage (warned partial result)",
    "no_assembly": "No assembly",
    "no_barcode": "No barcode recovered",
    "fail": "Failed coverage gate",
    "error": "Pipeline error",
}
STATUS_ORDER = list(STATUS_LABELS)

# Row/badge severity. `fail` and `error` share a colour deliberately —
# both are genuinely bad outcomes — but never share a label (spec
# §6a.4: "fail and error are not the same event... label them so an
# operator can tell them apart"), so the distinction survives in text
# even where it doesn't survive in colour. `low_coverage` gets its own
# `warning` colour, distinct from `ok`'s `success` (spec §6a.4's other
# requirement: a low_coverage row must not read as clean on a skim).
STATUS_SEVERITY = {
    "ok": "success",
    "low_coverage": "warning",
    "no_assembly": "danger",
    "no_barcode": "warning",
    "fail": "danger",
    "error": "danger",
}


def render(
    manifest: dict,
    template_dir: Path,
    static_dir: Path,
    out_path: Path,
) -> None:
    """Render a self-contained HTML report to `out_path`."""
    j2 = Environment(loader=FileSystemLoader(str(template_dir)))
    template = j2.get_template('run_index.html')

    context = build_context(manifest)
    context['static'] = get_static_file_contents(static_dir)

    rendered_html = template.render(**context)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered_html)


def build_context(manifest: dict) -> dict:
    breakdown = manifest.get('status_breakdown') or {}
    samplesheet = manifest.get('samplesheet') or {}
    reference_bundle = manifest.get('reference_bundle') or {}
    return {
        'title': 'WF5 run summary',
        'subtitle_html': REPORT_SUBTITLE_HTML,
        'sample_count': manifest.get('sample_count') or 0,
        'status_breakdown': [
            {
                'status': status,
                'label': STATUS_LABELS[status],
                'severity': STATUS_SEVERITY[status],
                'count': breakdown.get(status, 0),
            }
            for status in STATUS_ORDER
        ],
        'samples': [
            _sample_view(row) for row in manifest.get('samples') or []
        ],
        'pipeline_commit': manifest.get('pipeline_commit') or 'unknown',
        'invocation_timestamp': manifest.get('invocation_timestamp') or '-',
        'samplesheet_path': samplesheet.get('path') or '-',
        'samplesheet_sha256': samplesheet.get('sha256') or '-',
        'reference_bundle_version': reference_bundle.get('version'),
        'reference_bundle_generated_at': reference_bundle.get(
            'generated_at'),
    }


def _sample_view(row: dict) -> dict:
    status = row.get('sample_status')
    top_hit = row.get('top_blast_hit') or {}
    loci_total = row.get('n_barcodes_total')
    loci_passed = row.get('n_barcodes_passed')
    loci_partial = row.get('n_barcodes_partial') or 0
    coverage = row.get('coverage')
    return {
        **row,
        'status_label': STATUS_LABELS.get(status, status or '-'),
        'severity': STATUS_SEVERITY.get(status, 'secondary'),
        'coverage_text': f'{coverage}×' if coverage is not None else '-',
        'top_hit_text': (
            f"{top_hit['pident']:.2f}% — {top_hit.get('stitle')}"
            if top_hit.get('stitle') is not None
            and top_hit.get('pident') is not None
            else '-'
        ),
        'barcode_text': (
            (
                f'{loci_passed}/{loci_total}'
                + (f' (+{loci_partial} partial)' if loci_partial else '')
            ) if loci_total else '-'
        ),
        'report_link': f"{row.get('sample_id')}/report.html",
    }
