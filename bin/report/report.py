"""Render the per-sample `report.html` from `metadata.json`.

Single entry point, `render()`, called by `collate.py` (C6, task 42) —
wrapped there in a try/except so a rendering defect never takes the
sample's bundle down with it (task 43a §5.2).

`metadata.json` is the *only* JSON input (task 43a §3): sample_id,
kingdom, organelle, sample_status, coverage, assembly, homology,
barcodes, annotation and provenance all come from the parsed dict
passed in by the caller, never from globbing a result directory. The
handful of artefacts that are files rather than JSON fields — the
organelle map SVG, the Bandage graph PNG, the annotation GFF, the two
NanoPlot HTML reports — are passed in separately and inlined as
base64 `data:` URIs (or raw SVG markup) so the report stays a single
self-contained file (spec §6a.1).
"""

import os
from pathlib import Path
from typing import Optional

from jinja2 import Environment, FileSystemLoader

from .config import REPORT_SUBTITLE_HTML, REPORT_TITLE
from .filters.css_hash import css_hash
from .utils import file_data_uri, get_img_src

# Pipeline-order tab labels (spec §6a.2). Tab 1 (Overview) is rendered
# from its own template; 2-8 are registered here so the strip's shape
# is fixed once (task 43a §6) — 43b fills the panes.
STAGE_TABS = [
    "Sequencing data quality",
    "Recruitment & coverage gate",
    "De novo organelle assembly",
    "Assembly quality assessment",
    "Homology to reference databases",
    "Extracted barcode panel",
    "Annotated organelle genome map",
]

STATUS_LABELS = {
    "ok": "OK",
    "low_coverage": "Low coverage (warned partial result)",
    "no_assembly": "No assembly",
    "no_barcode": "No barcode recovered",
    "fail": "Failed coverage gate",
}


def render(
    metadata: dict,
    template_dir: Path,
    static_dir: Path,
    out_path: Path,
    params: Optional[dict] = None,
    nanoplot_raw_dir: Optional[Path] = None,
    nanoplot_clean_dir: Optional[Path] = None,
    organelle_map_svg: Optional[Path] = None,
    graph_png: Optional[Path] = None,
    annotation_gff: Optional[Path] = None,
) -> None:
    """Render a self-contained HTML report to `out_path`."""
    j2 = Environment(loader=FileSystemLoader(str(template_dir)))
    j2.filters['css_hash'] = css_hash
    template = j2.get_template('index.html')

    context = build_context(
        metadata,
        params=params or {},
        nanoplot_raw_dir=nanoplot_raw_dir,
        nanoplot_clean_dir=nanoplot_clean_dir,
        organelle_map_svg=organelle_map_svg,
        graph_png=graph_png,
        annotation_gff=annotation_gff,
    )
    context['static'] = get_static_file_contents(static_dir)

    rendered_html = template.render(**context)
    Path(out_path).write_text(rendered_html)


def get_static_file_contents(static_dir: Path) -> dict:
    """Return static file contents keyed for injection into the
    template — CSS/JS inlined as text, images as base64 `data:` URIs
    (spec §6a.1's self-containment requirement)."""
    static_files = {}
    for root, _, files in os.walk(static_dir):
        root = Path(root)
        if root.name == 'css':
            static_files['css'] = [
                f'/* {f} */\n' + (root / f).read_text()
                for f in sorted(files)
            ]
        elif root.name == 'js':
            static_files['js'] = [
                f'/* {f} */\n' + (root / f).read_text()
                for f in sorted(files)
            ]
        elif root.name == 'img':
            static_files['img'] = {
                f: get_img_src(root / f)
                for f in sorted(files)
            }
    return static_files


def build_context(
    metadata: dict,
    params: dict,
    nanoplot_raw_dir: Optional[Path],
    nanoplot_clean_dir: Optional[Path],
    organelle_map_svg: Optional[Path],
    graph_png: Optional[Path],
    annotation_gff: Optional[Path],
) -> dict:
    return {
        'title': REPORT_TITLE,
        'subtitle_html': REPORT_SUBTITLE_HTML,
        'sample_id': metadata.get('sample_id'),
        'metadata': metadata,
        'sample_status': metadata.get('sample_status'),
        'sample_status_label': STATUS_LABELS.get(
            metadata.get('sample_status'), metadata.get('sample_status')),
        'key_findings': key_findings(metadata),
        'warnings': build_warnings(metadata),
        'stage_tabs': STAGE_TABS,
        'parameters': params,
        'versions': (metadata.get('provenance') or {}).get(
            'tool_versions', {}),
        'organelle_map_svg': _read_svg(organelle_map_svg),
        'graph_png_src': _file_src(graph_png, 'image/png'),
        'annotation_gff_src': _file_src(annotation_gff, 'text/plain'),
        'nanoplot_raw_src': _nanoplot_report_src(nanoplot_raw_dir),
        'nanoplot_clean_src': _nanoplot_report_src(nanoplot_clean_dir),
    }


def _read_svg(path: Optional[Path]) -> Optional[str]:
    """The organelle map is inlined into the DOM, not an `<img>`, so
    43b's per-feature tooltips can bind to it. A missing or zero-byte
    SVG (today's stub — task 44 is out of scope here) renders as "not
    yet available" rather than raising (task 43a §4.3)."""
    if not path:
        return None
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return None
    return p.read_text()


def _file_src(path: Optional[Path], mime: str) -> Optional[str]:
    if not path:
        return None
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return None
    return file_data_uri(p, mime)


def _nanoplot_report_src(nanoplot_dir: Optional[Path]) -> Optional[str]:
    if not nanoplot_dir:
        return None
    path = Path(nanoplot_dir) / 'NanoPlot-report.html'
    if not path.is_file() or path.stat().st_size == 0:
        return None
    return file_data_uri(path, 'text/html')


# ---------------------------------------------------------------------------
# Key findings (spec §6a.1, CONSTITUTION.md principle 7) — the
# always-visible, negative-clarity-first summary. Every sample_status
# renders distinct text; no_assembly and no_barcode must never read as
# the same outcome, and low_coverage must never read as a failure
# (task 43a §7).
# ---------------------------------------------------------------------------


def key_findings(metadata: dict) -> list:
    status = metadata.get('sample_status')
    reason = metadata.get('sample_status_reason') or ''

    if status == 'fail':
        return [{
            'class': 'danger',
            'text': (
                'No organelle was assembled: the sample did not clear '
                f'the coverage gate’s hard minimum floor ({reason}).'
            ),
        }]
    if status == 'no_assembly':
        return [{
            'class': 'danger',
            'text': (
                'No organelle was assembled: recruited coverage cleared '
                'the gate, but no target contig was selected during '
                f'binning ({reason}).'
            ),
        }]

    findings = [
        _assembly_outcome_finding(metadata),
        _coverage_finding(metadata),
        _annotation_finding(metadata),
        _barcode_count_finding(metadata),
    ]
    top_hit = _top_blast_hit_finding(metadata)
    if top_hit:
        findings.append(top_hit)

    if status == 'no_barcode':
        findings.append({
            'class': 'warning',
            'text': (
                'An organelle genome was assembled and annotated, but no '
                f'barcode locus passed validation ({reason}).'
            ),
        })
    elif status == 'low_coverage':
        findings.append({
            'class': 'info',
            'text': (
                'This is a real recovery produced below the coverage '
                'warn floor — read it as a partial result with the '
                'caveats in the Recruitment & coverage gate tab, not as '
                'a degraded or failed one.'
            ),
        })
    else:
        findings.append({'class': 'success', 'text': 'Clean recovery.'})
    return findings


def _assembly_outcome_finding(metadata: dict) -> dict:
    assembly = metadata.get('assembly') or {}
    n = assembly.get('contig_count')
    total_bp = assembly.get('total_bp')
    if n and total_bp:
        text = f'Assembled {n} contig(s) totalling {total_bp:,} bp.'
    else:
        text = 'An organelle assembly was produced.'
    return {'class': 'success', 'text': text}


def _coverage_finding(metadata: dict) -> dict:
    gate = (metadata.get('coverage') or {}).get('gate') or {}
    cov = gate.get('estimated_cov')
    if cov is None:
        text = 'Coverage estimate unavailable.'
    else:
        text = f'Estimated recruited coverage: {cov}×.'
    return {'class': 'secondary', 'text': text}


def _annotation_finding(metadata: dict) -> dict:
    counts = (metadata.get('annotation') or {}).get('feature_counts') or {}
    gene_n = counts.get('gene')
    if gene_n is None:
        text = 'No annotation was produced.'
    else:
        text = f'{gene_n} gene feature(s) annotated.'
    return {'class': 'secondary', 'text': text}


def _barcode_count_finding(metadata: dict) -> dict:
    barcodes = metadata.get('barcodes') or {}
    loci = barcodes.get('loci') or []
    n_passed = barcodes.get('n_passed') or 0
    if not loci:
        text = 'No barcode loci were evaluated.'
    else:
        text = f'{n_passed}/{len(loci)} panel loci passed barcode validation.'
    return {'class': 'secondary', 'text': text}


def _top_blast_hit_finding(metadata: dict) -> Optional[dict]:
    hits = (metadata.get('homology') or {}).get('top_hits') or []
    if not hits:
        return None
    best = max(hits, key=lambda h: h.get('bitscore', 0))
    return {
        'class': 'secondary',
        'text': (
            f"Top BLAST hit: {best.get('stitle')} "
            f"({best.get('pident')}% identity)."
        ),
    }


# ---------------------------------------------------------------------------
# Warnings mirror (task 43a §6) — every entry carries text, severity
# and the tab it belongs to. Overview iterates this list so a warning
# never exists only behind the tab that owns it; 43b appends to it
# rather than inventing per-tab warning markup.
# ---------------------------------------------------------------------------


def build_warnings(metadata: dict) -> list:
    warnings = []
    status = metadata.get('sample_status')

    if status == 'low_coverage':
        warnings.append({
            'severity': 'warning',
            'tab': 'Recruitment & coverage gate',
            'text': (
                'This sample assembled below the coverage warn floor. '
                'Missing genes and fragmented contigs are expected, not '
                'a sign that something went wrong.'
            ),
        })

    gate = (metadata.get('coverage') or {}).get('gate') or {}
    if gate.get('coverage_basis') == 'total_recruited':
        warnings.append({
            'severity': 'info',
            'tab': 'Recruitment & coverage gate',
            'text': (
                'The sibling-organelle split was unavailable for this '
                'sample — the coverage estimate may over-state '
                'target depth (spec §2.1.6).'
            ),
        })

    annotation = metadata.get('annotation') or {}
    if annotation.get('status') == 'annotator_failed':
        warnings.append({
            'severity': 'danger',
            'tab': 'Annotated organelle genome map',
            'text': (
                f"The non-CDS annotator failed to run: "
                f"{annotation.get('reason')}"
            ),
        })
    if annotation.get('genetic_code_agreement') is False:
        warnings.append({
            'severity': 'warning',
            'tab': 'Annotated organelle genome map',
            'text': (
                'ANNOTATE and EXTRACT_BARCODES selected different '
                'genetic-code tables for this sample — see the '
                'provenance panel.'
            ),
        })
    crosscheck = annotation.get('cds_crosscheck') or {}
    if crosscheck.get('annotator_only'):
        warnings.append({
            'severity': 'warning',
            'tab': 'Annotated organelle genome map',
            'text': (
                'miniprot and the non-CDS annotator disagree on one or '
                'more genes — see the cross-check table.'
            ),
        })

    if status == 'no_barcode':
        warnings.append({
            'severity': 'warning',
            'tab': 'Extracted barcode panel',
            'text': 'No barcode locus passed validation for this sample.',
        })

    return warnings
