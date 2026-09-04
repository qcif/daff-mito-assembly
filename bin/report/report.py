"""Render the per-sample `report.html` from `metadata.json`.

Single entry point, `render()`, called by `collate.py` (C6, task 42) —
wrapped there in a try/except so a rendering defect never takes the
sample's bundle down with it (task 43a §5.2).

`metadata.json` is the *only* JSON input (task 43a §3): sample_id,
kingdom, organelle, sample_status, coverage, assembly, homology,
barcodes, annotation and provenance all come from the parsed dict
passed in by the caller, never from globbing a result directory. The
handful of artefacts that are files rather than JSON fields — the
organelle map SVG, the Bandage graph PNG, the annotation GFF, the
recovered-barcode FASTA, the two NanoPlot HTML reports — are passed in
separately and inlined as base64 `data:` URIs (or raw SVG markup) so
the report stays a single self-contained file (spec §6a.1).

Four reader-perspective tabs (task 43b, spec §6a.2; reordered task 46
§4), not one tab per pipeline stage: Overview, Validation, Assembly,
Barcodes. Overview is the default and mirrors every warning raised on
any other tab, via a `warnings` list in the render context, so a
caveat is never visible only behind a tab a reader might not open.
"""

import os
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

from jinja2 import Environment, FileSystemLoader

from . import confidence
from .config import REPORT_SUBTITLE_HTML
from .filters.css_hash import css_hash
from .utils import file_data_uri, get_img_src

# Organelle-specific H1 title (task 43b §2.1 item 1, spec §6a.1) — a
# constant here, not template logic, so "which organelle am I looking
# at" never has to be inferred from the body.
ORGANELLE_TITLES = {
    "mt": "Mitochondrial genome assembly",
    "pt": "Chloroplast genome assembly",
}
DEFAULT_TITLE = "Organelle genome assembly"

# Full organelle names (task 46 §5 item 3) — the single source of
# truth for "which organelle" wherever the report spells it out in
# prose, not just the Input params row it started from. "Mitochondrion"
# is the singular (one organelle genome per report), confirmed with
# the project owner over the review's plural "Mitochondria".
ORGANELLE_NAMES = {
    "mt": "Mitochondrion",
    "pt": "Chloroplast",
}

# Four reader-perspective tabs (spec §6a.2, task 46 §4) — the
# pipeline's stage boundaries are an implementation detail; a
# biosecurity officer asks "did it work / should I trust it / what did
# I get / which barcodes can I use", in that order. Validation precedes
# Assembly: trustworthiness precedes content, and on a soft-failed
# sample Validation is already the terminal tab (43b §3.3), so it sits
# where the reader arrives next after Overview.
TAB_NAMES = ["Overview", "Validation", "Assembly", "Barcodes"]

# One ordered list of pane templates, driven by `TAB_NAMES` itself
# (task 46 §4), so the tab-button strip and the pane markup cannot
# drift out of sync the way the previous hardcoded `tab-1`/`tab-2`
# includes did.
TAB_TEMPLATES = {
    "Overview": "components/overview.html",
    "Validation": "components/validation.html",
    "Assembly": "components/assembly.html",
    "Barcodes": "components/barcodes.html",
}

STATUS_LABELS = {
    "ok": "OK",
    "low_coverage": "Low coverage (warned partial result)",
    "no_assembly": "No assembly",
    "no_barcode": "No barcode recovered",
    "fail": "Failed coverage gate",
}

# A sample that never reached assembly — the Assembly and Barcodes
# tabs must say so rather than render empty tables (spec §6a.2, task
# 43b §3.3).
TERMINAL_STATUSES = ("fail", "no_assembly")

# Species-level identity floor for the Validation tab's below-threshold
# flag (task 43b §5.2) — a commonly used COI/barcoding species-level
# convention. Provisional pending spec §9's benchmarking sweep; the
# tab only ever *flags* against it, never *assigns* a taxon.
SPECIES_IDENTITY_THRESHOLD = 97.0

BARCODE_DROPOUT_REASONS = {
    "not_found": "Locus not found in the annotated CDS features.",
    "invalid_length": "Extracted sequence length is not a multiple of 3.",
    "identity_below_floor": (
        "Protein identity to the miniprot reference is below the "
        "validation floor."
    ),
    "internal_stop_codon": (
        "The extracted ORF carries an internal stop codon under every "
        "configured genetic-code table."
    ),
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
    barcodes_fasta: Optional[Path] = None,
    workflow_start: Optional[str] = None,
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
        barcodes_fasta=barcodes_fasta,
        workflow_start=workflow_start,
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
    barcodes_fasta: Optional[Path],
    workflow_start: Optional[str],
) -> dict:
    status = metadata.get('sample_status')
    title = ORGANELLE_TITLES.get(metadata.get('organelle'), DEFAULT_TITLE)
    return {
        'title': title,
        'subtitle_html': REPORT_SUBTITLE_HTML,
        'sample_id': metadata.get('sample_id'),
        'metadata': metadata,
        'organelle_name': ORGANELLE_NAMES.get(
            metadata.get('organelle'), metadata.get('organelle') or '-'),
        'sample_status': status,
        'sample_status_label': STATUS_LABELS.get(status, status),
        'terminal': status in TERMINAL_STATUSES,
        'key_findings': key_findings(metadata),
        'warnings': build_warnings(metadata),
        'tab_names': TAB_NAMES,
        'tab_templates': [TAB_TEMPLATES[name] for name in TAB_NAMES],
        'parameters': params,
        'facility': params.get('facility') or '-',
        'analyst_name': params.get('analyst_name') or '-',
        'wall_time': wall_time_context(workflow_start),
        'versions': (metadata.get('provenance') or {}).get(
            'tool_versions', {}),
        'organelle_map_svg': _read_svg(organelle_map_svg),
        'graph_png_src': _file_src(graph_png, 'image/png'),
        'annotation_gff_src': _file_src(annotation_gff, 'text/plain'),
        'barcodes_fasta_src': _file_src(barcodes_fasta, 'text/plain'),
        'nanoplot_raw_src': _nanoplot_report_src(nanoplot_raw_dir),
        'nanoplot_clean_src': _nanoplot_report_src(nanoplot_clean_dir),
        'assembly_view': assembly_view(metadata),
        'validation_view': validation_view(metadata),
        'barcodes_view': barcodes_view(metadata, barcodes_fasta),
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
# Wall time (task 43b §2.1 item 3) — restored from the deleted
# boilerplate component with a real data source: `workflow.start`
# passed through COLLATE as a CLI argument, end time taken at render.
# ---------------------------------------------------------------------------


def wall_time_context(workflow_start: Optional[str]) -> dict:
    if not workflow_start:
        return {'start_time': None, 'end_time': None, 'duration': None}
    try:
        start = datetime.fromisoformat(workflow_start)
    except ValueError:
        return {'start_time': workflow_start, 'end_time': None,
                'duration': None}
    end = datetime.now()
    delta = max(end - start, timedelta(0))
    total_seconds = int(delta.total_seconds())
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return {
        'start_time': start.strftime('%Y-%m-%d %H:%M:%S'),
        'end_time': end.strftime('%Y-%m-%d %H:%M:%S'),
        'duration': f'{hours:02d}:{minutes:02d}:{seconds:02d}',
    }


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
            'label': 'Outcome',
            'text': (
                'No organelle was assembled: the sample did not clear '
                f'the coverage gate’s hard minimum floor ({reason}).'
            ),
        }]
    if status == 'no_assembly':
        return [{
            'class': 'danger',
            'label': 'Outcome',
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
            'label': 'Outcome',
            'text': (
                'An organelle genome was assembled and annotated, but no '
                f'barcode locus passed validation ({reason}).'
            ),
        })
    elif status == 'low_coverage':
        findings.append({
            'class': 'info',
            'label': 'Outcome',
            'text': (
                'This is a real recovery produced below the coverage '
                'warn floor — read it as a partial result with the '
                'caveats on the Validation tab, not as a degraded or '
                'failed one.'
            ),
        })
    else:
        findings.append({
            'class': 'success', 'label': 'Outcome', 'text': 'Clean recovery.',
        })
    return findings


def _assembly_outcome_finding(metadata: dict) -> dict:
    assembly = metadata.get('assembly') or {}
    n = assembly.get('contig_count')
    total_bp = assembly.get('total_bp')
    if n and total_bp:
        text = f'Assembled {n} contig(s) totalling {total_bp:,} bp.'
    else:
        text = 'An organelle assembly was produced.'
    return {'class': 'success', 'label': 'Assembly', 'text': text}


def _coverage_finding(metadata: dict) -> dict:
    gate = (metadata.get('coverage') or {}).get('gate') or {}
    cov = gate.get('estimated_cov')
    if cov is None:
        text = 'Coverage estimate unavailable.'
    else:
        text = f'Estimated recruited coverage: {cov}×.'
    return {'class': 'secondary', 'label': 'Coverage', 'text': text}


def _annotation_finding(metadata: dict) -> dict:
    counts = (metadata.get('annotation') or {}).get('feature_counts') or {}
    gene_n = counts.get('gene')
    if gene_n is None:
        text = 'No annotation was produced.'
    else:
        text = f'{gene_n} gene feature(s) annotated.'
    return {'class': 'secondary', 'label': 'Annotation', 'text': text}


def _barcode_count_finding(metadata: dict) -> dict:
    barcodes = metadata.get('barcodes') or {}
    loci = barcodes.get('loci') or []
    n_passed = barcodes.get('n_passed') or 0
    if not loci:
        text = 'No barcode loci were evaluated.'
    else:
        text = f'{n_passed}/{len(loci)} panel loci passed barcode validation.'
    return {'class': 'secondary', 'label': 'Barcodes', 'text': text}


def _top_blast_hit_finding(metadata: dict) -> Optional[dict]:
    hits = (metadata.get('homology') or {}).get('top_hits') or []
    if not hits:
        return None
    best = max(hits, key=lambda h: h.get('bitscore', 0))
    return {
        'class': 'secondary',
        'label': 'Top hit',
        'text': (
            f"{best.get('stitle')} ({best.get('pident')}% identity)."
        ),
    }


# ---------------------------------------------------------------------------
# Warnings mirror (task 43a §6, extended task 43b §3) — every entry
# carries text, severity and the tab it belongs to. Overview iterates
# this list so a warning never exists only behind the tab that owns
# it.
# ---------------------------------------------------------------------------


def build_warnings(metadata: dict) -> list:
    warnings = []
    status = metadata.get('sample_status')
    kingdom = metadata.get('kingdom')

    if status == 'low_coverage':
        warnings.append({
            'severity': 'warning',
            'tab': 'Validation',
            'text': (
                'This sample assembled below the coverage warn floor. '
                'Missing genes and fragmented contigs are expected, not '
                'a sign that something went wrong.'
            ),
        })

    gate = (metadata.get('coverage') or {}).get('gate') or {}
    if gate.get('coverage_basis') == 'total_recruited':
        # Task 46 §8 item 3 — 'warning', not 'info': this caveats a
        # number that is otherwise reported at face value, and the
        # Overview mirror must agree with the Validation tab it
        # mirrors (rule 18).
        warnings.append({
            'severity': 'warning',
            'tab': 'Validation',
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
            'tab': 'Assembly',
            'text': (
                f"The non-CDS annotator failed to run: "
                f"{annotation.get('reason')}"
            ),
        })
    if annotation.get('genetic_code_agreement') is False:
        warnings.append({
            'severity': 'warning',
            'tab': 'Validation',
            'text': (
                'ANNOTATE and EXTRACT_BARCODES selected different '
                'genetic-code tables for this sample — see the '
                'annotation cross-check section.'
            ),
        })
    crosscheck = annotation.get('cds_crosscheck') or {}
    if crosscheck.get('annotator_only'):
        warnings.append({
            'severity': 'warning',
            'tab': 'Validation',
            'text': (
                'miniprot and the non-CDS annotator disagree on one or '
                'more genes — see the cross-check table.'
            ),
        })

    if kingdom == 'plant' and annotation.get('status') == 'ok_cds_only':
        warnings.append({
            'severity': 'info',
            'tab': 'Assembly',
            'text': (
                'This annotation is CDS-only: tRNA/rRNA calling is not '
                'yet available for plant targets (spec §8 item 3).'
            ),
        })

    plastid = (
        (metadata.get('assembly') or {}).get('bin_metadata') or {}
    ).get('plastid_canonicalisation') or {}
    if plastid.get('branch') == 'canonical' and not plastid.get(
            'substitution_applied'):
        warnings.append({
            'severity': 'warning',
            'tab': 'Assembly',
            'text': (
                'The assembly graph shows a canonical quadripartite '
                'structure, but no taxonomic support — the substitution '
                'was withheld: '
                f"{plastid.get('substitution_withheld_reason')}"
            ),
        })

    if status == 'no_barcode':
        warnings.append({
            'severity': 'warning',
            'tab': 'Barcodes',
            'text': 'No barcode locus passed validation for this sample.',
        })

    top_hits = (metadata.get('homology') or {}).get('top_hits') or []
    if any(
        h.get('pident') is not None
        and h['pident'] < SPECIES_IDENTITY_THRESHOLD
        for h in top_hits
    ):
        warnings.append({
            'severity': 'info',
            'tab': 'Validation',
            'text': (
                'Top-hit identity is below the species-level threshold '
                f'({SPECIES_IDENTITY_THRESHOLD}%) — this flags a novel '
                'or under-represented taxon; it does not assign one.'
            ),
        })
    return warnings


# ---------------------------------------------------------------------------
# Assembly tab (task 43b §3.2)
# ---------------------------------------------------------------------------


def assembly_view(metadata: dict) -> dict:
    assembly = metadata.get('assembly') or {}
    bin_metadata = assembly.get('bin_metadata') or {}
    classifications = {
        c.get('contig_id'): c.get('classification')
        for c in bin_metadata.get('contigs') or []
    }
    contigs = []
    for c in assembly.get('contigs') or []:
        cls = classifications.get(c.get('contig'))
        contigs.append({
            **c,
            'bucket': _contig_bucket(cls),
        })

    annotation = metadata.get('annotation') or {}
    cds_scores = sorted(
        annotation.get('cds_scores') or [],
        key=lambda s: (s.get('seqid') or '', s.get('start') or 0),
    )
    scored = [_score_row(s) for s in cds_scores]

    return {
        'contigs': contigs,
        'coverage_chart': _coverage_chart_data(contigs),
        'plastid': (bin_metadata or {}).get('plastid_canonicalisation'),
        'target_source': bin_metadata.get('target_source'),
        'cds_scores': scored,
        'cds_only': (
            metadata.get('kingdom') == 'plant'
            and annotation.get('status') == 'ok_cds_only'
        ),
        'annotator_failed': annotation.get('status') == 'annotator_failed',
    }


def _contig_bucket(classification: Optional[str]) -> str:
    if classification == 'target_candidate':
        return 'target'
    if classification == 'secondary_target':
        return 'secondary'
    return 'off-target'


_BUCKET_COLOURS = {
    'target': '#2ca02c',
    'secondary': '#ff7f0e',
    'off-target': '#7f7f7f',
}


def _coverage_chart_data(contigs: list) -> dict:
    contigs = [c for c in contigs if c.get('coverage') is not None]
    return {
        'x': [c.get('contig') for c in contigs],
        'y': [c.get('coverage') for c in contigs],
        'colors': [_BUCKET_COLOURS[c['bucket']] for c in contigs],
    }


def _score_row(score: dict) -> dict:
    quadrant = confidence.classify(score.get('pident'), score.get('qcovhsp'))
    return {
        **score,
        'quadrant': quadrant,
        'quadrant_label': confidence.label(quadrant),
        'quadrant_description': confidence.description(quadrant),
        'completeness_severity': confidence.axis_severity(
            quadrant, confidence.COMPLETENESS),
        'completeness_tooltip': confidence.axis_tooltip(
            quadrant, confidence.COMPLETENESS),
        'reference_severity': confidence.axis_severity(
            quadrant, confidence.REFERENCE),
        'reference_tooltip': confidence.axis_tooltip(
            quadrant, confidence.REFERENCE),
    }


# ---------------------------------------------------------------------------
# Validation tab (task 43b §3.3)
# ---------------------------------------------------------------------------


def validation_view(metadata: dict) -> dict:
    coverage = metadata.get('coverage') or {}
    gate = coverage.get('gate') or {}
    estimate = coverage.get('estimate') or {}
    homology = metadata.get('homology') or {}
    annotation = metadata.get('annotation') or {}
    top_hits = [
        {**h, 'below_species_threshold': (
            h.get('pident') is not None
            and h['pident'] < SPECIES_IDENTITY_THRESHOLD
        )}
        for h in homology.get('top_hits') or []
    ]
    return {
        'gate': gate,
        'estimate': estimate,
        'recruitment': coverage.get('recruitment'),
        'top_hits': top_hits,
        'cds_crosscheck': annotation.get('cds_crosscheck'),
        'genetic_code_annotate': annotation.get('genetic_code_annotate'),
        'genetic_code_cds': annotation.get('genetic_code_cds'),
        'genetic_code_agreement': annotation.get('genetic_code_agreement'),
        'flye_depth': _flye_depth(metadata),
    }


def _flye_depth(metadata: dict) -> Optional[dict]:
    """Flye's own mean coverage on the contig(s) actually selected as
    the target (task 46 §9) — read against `bin_metadata.
    contigs_selected`, distinct from the gate's recruited-read
    estimate above it. `None` for `fail` / `no_assembly`, which have
    no selected contig — the template renders that as `-`, never `0`
    (a measured zero depth would be a false claim)."""
    bin_metadata = (metadata.get('assembly') or {}).get(
        'bin_metadata') or {}
    selected = bin_metadata.get('contigs_selected') or []
    if not selected:
        return None
    coverage_by_contig = {
        c.get('contig'): c.get('coverage')
        for c in (metadata.get('assembly') or {}).get('contigs') or []
    }
    return {
        'contigs': selected,
        'coverages': [coverage_by_contig.get(cid) for cid in selected],
    }


# ---------------------------------------------------------------------------
# Barcodes tab (task 43b §3.4)
# ---------------------------------------------------------------------------


def barcodes_view(
    metadata: dict, barcodes_fasta: Optional[Path] = None,
) -> dict:
    barcodes = metadata.get('barcodes') or {}
    loci = barcodes.get('loci') or []
    sequences = _parse_fasta(barcodes_fasta)
    passed = [
        {**locus, 'sequence': sequences.get(_barcode_id(locus), '')}
        for locus in loci if locus.get('status') == 'pass'
    ]
    dropped = [
        {
            **locus,
            'reason_text': BARCODE_DROPOUT_REASONS.get(
                locus.get('reason'), locus.get('reason') or 'Unknown.'),
        }
        for locus in loci if locus.get('status') != 'pass'
    ]
    return {
        'passed': passed,
        'dropped': dropped,
        'n_passed': barcodes.get('n_passed'),
        'n_total': len(loci),
    }


def _barcode_id(locus: dict) -> str:
    return (
        f"{locus.get('gene')}_{locus.get('seqid')}_"
        f"{locus.get('start')}_{locus.get('end')}"
    )


def _parse_fasta(path: Optional[Path]) -> dict:
    if not path:
        return {}
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return {}
    sequences = {}
    current_id = None
    chunks: list = []
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            if current_id is not None:
                sequences[current_id] = ''.join(chunks)
            current_id = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if current_id is not None:
        sequences[current_id] = ''.join(chunks)
    return sequences
