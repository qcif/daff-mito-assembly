#!/usr/bin/env python3
"""Stage 14 -- annotated organelle map (task 44).

In-house Jinja+SVG renderer, chosen over pyCirclize/plotgenes/Circos
(see task 44_organelle_map.md §1): spec §6a.1 wants clickable/hoverable
gene features carrying tooltip metadata, which means an inline,
DOM-injectable SVG rather than a rendered picture -- a circos-style
library would still need a tooltip layer built on top of it, and the
in-house route reuses the C6/C7 report container instead of adding a
fourth image (spec §6a.5).

Reads the merged annotation GFF3 + ``annotation_summary.json`` (task 40)
and ``bin_metadata.json`` (C3/C4) and lays the annotated features out on
a circular genome diagram. Tooltips use a native SVG ``<title>`` child
per feature -- no JavaScript, consistent with the report's
self-contained-file requirement (spec §6a.5, and the same choice task
47 makes for the assembly graph).

A single physical gene locus routinely surfaces as several near-
identical raw GFF records, because the broad miniprot pass (spec §8
item 3) independently hits every reference-panel ortholog at that
locus. Genomically-overlapping same-strand CDS calls are collapsed onto
one drawn arc, keeping the best-scored call as the primary annotation
and listing any other gene labels seen at that locus in the tooltip
(``alt_genes``) rather than silently dropping them (CONSTITUTION
principle 7). tRNA/rRNA calls have only one source (MITOS2) and are
never clustered.

``ORGANELLE_MAP`` is the sole stage aware of the plastid quadripartite
structure (spec §3.6 step 5): when ``bin_metadata.json`` records a
canonical plastid substitution, a second panel is rendered for
``path2`` by reflecting every feature that falls inside the SSC segment
(the segment whose orientation differs between the two isoforms) and
flipping its strand; a feature that spans the SSC boundary cannot be
cleanly remapped and is dropped from the path2 panel only.

Never raises past ``main()`` -- a rendering defect must not fail the
sample (CONSTITUTION principle 8); on any error the process still exits
0, leaving an empty output file for the report's existing "not yet
available" fallback (``bin/report/report.py``'s ``_read_svg``).
"""

import argparse
import json
import math
import sys
from pathlib import Path
from xml.sax.saxutils import escape

from annotation_gff import parse_annotation_gff

RING_RADIUS = 200
FEATURE_RING_WIDTH = 20
FEATURE_RING_GAP = 6
PANEL_MARGIN = 70
MIN_WEDGE_RADIANS = 0.01
LEGEND_ROW_HEIGHT = 18

KIND_COLOURS = {
    'CDS': '#4C78A8',
    'tRNA': '#F58518',
    'rRNA': '#54A24B',
}
DEFAULT_COLOUR = '#999999'
BACKBONE_COLOUR = '#cccccc'


# ── input parsing ────────────────────────────────────────────────────

def parse_non_cds_features(gff_path) -> list:
    """tRNA / rRNA top-level rows. ``annotation_gff.parse_annotation_cds``
    only covers the mRNA/CDS hierarchy (task 40's shared parser), so
    this stage reads the remaining feature types itself."""
    features = []
    with open(gff_path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) != 9:
                continue
            seqid, source, ftype, start, end, _score, strand, \
                _phase, attr = fields
            if ftype not in ('tRNA', 'rRNA'):
                continue
            attrs = {}
            for kv in attr.strip().split(';'):
                if kv and '=' in kv:
                    key, value = kv.split('=', 1)
                    attrs[key] = value
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            features.append({
                'seqid': seqid, 'source': source, 'kind': ftype,
                'strand': strand, 'start': start_i, 'end': end_i,
                'gene': attrs.get('Name') or attrs.get('gene_id') or ftype,
                'pident': None, 'qcovhsp': None, 'bitscore': None,
            })
    return features


def load_cds_scores(annotation_summary_path) -> dict:
    if not annotation_summary_path:
        return {}
    path = Path(annotation_summary_path)
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    try:
        summary = json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}
    return {
        row['id']: row for row in summary.get('cds_scores', [])
        if row.get('id')
    }


def cds_features(gff_path, cds_scores: dict) -> list:
    features = []
    for rec in parse_annotation_gff(gff_path):
        score = cds_scores.get(rec['id'], {})
        features.append({
            'seqid': rec['seqid'], 'source': rec['source'], 'kind': 'CDS',
            'strand': rec['strand'], 'start': rec['start'],
            'end': rec['end'], 'gene': score.get('gene', rec['gene']),
            'pident': score.get('pident'), 'qcovhsp': score.get('qcovhsp'),
            'bitscore': score.get('bitscore'),
        })
    return features


def load_bin_metadata(bin_metadata_path) -> dict:
    if not bin_metadata_path:
        return {}
    path = Path(bin_metadata_path)
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}


# ── locus clustering ─────────────────────────────────────────────────

def _score_key(feature: dict):
    return (
        feature.get('pident') is not None,
        feature.get('bitscore') or 0,
        feature['end'] - feature['start'],
    )


def _cluster_group(group: list) -> list:
    group = sorted(group, key=lambda f: f['start'])
    clusters, current, current_end = [], [], None
    for feature in group:
        if current and feature['start'] <= current_end:
            current.append(feature)
            current_end = max(current_end, feature['end'])
        else:
            if current:
                clusters.append(current)
            current, current_end = [feature], feature['end']
    if current:
        clusters.append(current)
    return clusters


def cluster_loci(cds_feats: list, other_feats: list) -> list:
    """One drawn locus per singleton non-CDS feature, and one per
    genomically-overlapping same-(seqid, strand) cluster of CDS
    features -- see module docstring."""
    loci = [{**f, 'alt_genes': []} for f in other_feats]

    by_seqid_strand = {}
    for feature in cds_feats:
        by_seqid_strand.setdefault(
            (feature['seqid'], feature['strand']), []).append(feature)

    for group in by_seqid_strand.values():
        for cluster in _cluster_group(group):
            primary = max(cluster, key=_score_key)
            alt_genes = sorted({
                f['gene'] for f in cluster if f['gene'] != primary['gene']
            })
            loci.append({**primary, 'alt_genes': alt_genes})

    return loci


# ── genome length + plastid path2 remap ─────────────────────────────

def genome_length(metadata: dict, loci: list):
    canon = metadata.get('plastid_canonicalisation') or {}
    if canon.get('substitution_applied') and canon.get('path1_len'):
        return canon['path1_len']

    contigs_selected = metadata.get('contigs_selected') or []
    contigs = {
        c['contig_id']: c for c in metadata.get('contigs', [])
        if c.get('contig_id')
    }
    if contigs_selected and contigs_selected[0] in contigs:
        length = contigs[contigs_selected[0]].get('length_bp')
        if length:
            return length

    if loci:
        return max(locus['end'] for locus in loci)
    return None


def plastid_path2_segments(metadata: dict):
    """Returns (lsc_len, ir_len, ssc_len) if ``bin_metadata.json``
    records a canonical plastid substitution with segment lengths,
    else ``None`` -- the trigger for rendering a second, path2 panel
    (spec §3.6 step 5)."""
    canon = metadata.get('plastid_canonicalisation') or {}
    if not canon.get('substitution_applied'):
        return None
    lsc_len, ir_len, ssc_len = (
        canon.get('lsc_len'), canon.get('ir_len'), canon.get('ssc_len'))
    if not all(isinstance(v, int) for v in (lsc_len, ir_len, ssc_len)):
        return None
    return lsc_len, ir_len, ssc_len


_STRAND_FLIP = {'+': '-', '-': '+'}


def remap_to_path2(loci: list, segments: tuple) -> list:
    """path2 = LSC + IR + revcomp(SSC) + revcomp(IR) -- identical to
    path1 outside the SSC segment, whose bases (and hence any feature
    inside it) are reversed in place (spec/plastid-canonicalisation.md).
    A feature straddling the SSC boundary has no valid position in
    path2's coordinate frame and is dropped from this panel only."""
    lsc_len, ir_len, ssc_len = segments
    ssc_start = lsc_len + ir_len  # 0-based bp preceding the SSC segment
    ssc_end = ssc_start + ssc_len  # last 1-based SSC position

    remapped = []
    for locus in loci:
        start, end = locus['start'], locus['end']
        if end <= ssc_start or start > ssc_end:
            remapped.append(locus)
            continue
        if start <= ssc_start or end > ssc_end:
            continue
        new_start = 2 * ssc_start + ssc_len + 1 - end
        new_end = 2 * ssc_start + ssc_len + 1 - start
        remapped.append({
            **locus, 'start': new_start, 'end': new_end,
            'strand': _STRAND_FLIP.get(locus['strand'], locus['strand']),
        })
    return remapped


# ── SVG rendering ─────────────────────────────────────────────────────

def _point(cx, cy, r, theta):
    return cx + r * math.cos(theta), cy + r * math.sin(theta)


def _angle(pos, length):
    return -math.pi / 2 + 2 * math.pi * (pos / length)


def _wedge_path(cx, cy, r_inner, r_outer, theta1, theta2) -> str:
    if theta2 - theta1 < MIN_WEDGE_RADIANS:
        theta2 = theta1 + MIN_WEDGE_RADIANS
    large_arc = 1 if (theta2 - theta1) > math.pi else 0
    x1o, y1o = _point(cx, cy, r_outer, theta1)
    x2o, y2o = _point(cx, cy, r_outer, theta2)
    x2i, y2i = _point(cx, cy, r_inner, theta2)
    x1i, y1i = _point(cx, cy, r_inner, theta1)
    return (
        f'M {x1o:.2f} {y1o:.2f} '
        f'A {r_outer:.2f} {r_outer:.2f} 0 {large_arc} 1 {x2o:.2f} {y2o:.2f} '
        f'L {x2i:.2f} {y2i:.2f} '
        f'A {r_inner:.2f} {r_inner:.2f} 0 {large_arc} 0 {x1i:.2f} {y1i:.2f} '
        f'Z'
    )


def _ring_radii(strand: str) -> tuple:
    if strand == '-':
        return RING_RADIUS - FEATURE_RING_GAP - FEATURE_RING_WIDTH, \
            RING_RADIUS - FEATURE_RING_GAP
    return RING_RADIUS + FEATURE_RING_GAP, \
        RING_RADIUS + FEATURE_RING_GAP + FEATURE_RING_WIDTH


def _tooltip_text(locus: dict) -> str:
    lines = [
        locus['gene'],
        f"{locus['kind']} | {locus['seqid']}:{locus['start']}-"
        f"{locus['end']} ({locus['strand']})",
        f"source: {locus['source']}",
    ]
    if locus.get('alt_genes'):
        lines.append('also called: ' + ', '.join(locus['alt_genes']))
    if locus.get('pident') is not None:
        lines.append(
            f"pident={locus['pident']} qcovhsp={locus['qcovhsp']} "
            f"bitscore={locus['bitscore']}"
        )
    return '\n'.join(lines)


def _render_feature(cx, cy, locus: dict, length: int) -> str:
    theta1 = _angle(locus['start'] - 1, length)
    theta2 = _angle(locus['end'], length)
    r_inner, r_outer = _ring_radii(locus['strand'])
    colour = KIND_COLOURS.get(locus['kind'], DEFAULT_COLOUR)
    path_d = _wedge_path(cx, cy, r_inner, r_outer, theta1, theta2)
    gene = escape(locus['gene'])
    tooltip = escape(_tooltip_text(locus))
    return (
        f'<g data-gene="{gene}" data-kind="{locus["kind"]}" '
        f'data-source="{escape(locus["source"])}">'
        f'<title>{tooltip}</title>'
        f'<path d="{path_d}" fill="{colour}" stroke="#ffffff" '
        f'stroke-width="0.5"/></g>'
    )


def _render_panel(cx, cy, loci: list, length: int, label: str) -> str:
    parts = [
        f'<circle cx="{cx}" cy="{cy}" r="{RING_RADIUS}" fill="none" '
        f'stroke="{BACKBONE_COLOUR}" stroke-width="1.5"/>',
        f'<text x="{cx}" y="{cy}" text-anchor="middle" '
        f'dominant-baseline="middle" font-size="13">{escape(label)}'
        f'<tspan x="{cx}" dy="16">{length:,} bp</tspan></text>',
    ]
    for locus in loci:
        parts.append(_render_feature(cx, cy, locus, length))
    return ''.join(parts)


def _render_legend(x, y) -> str:
    rows = []
    for i, (kind, colour) in enumerate(KIND_COLOURS.items()):
        row_y = y + i * LEGEND_ROW_HEIGHT
        rows.append(
            f'<rect x="{x}" y="{row_y}" width="12" height="12" '
            f'fill="{colour}"/>'
            f'<text x="{x + 18}" y="{row_y + 10}" font-size="12">'
            f'{kind}</text>'
        )
    rows.append(
        f'<text x="{x}" y="{y + len(KIND_COLOURS) * LEGEND_ROW_HEIGHT + 14}"'
        f' font-size="11" fill="#666666">outer ring: + strand / inner '
        f'ring: - strand</text>'
    )
    return ''.join(rows)


def render_svg(panels: list) -> str:
    """``panels`` is a list of ``(label, loci, length)`` -- one entry
    for the primary map, two when the plastid path2 panel applies."""
    panel_width = 2 * (RING_RADIUS + PANEL_MARGIN)
    height = 2 * (RING_RADIUS + PANEL_MARGIN) + LEGEND_ROW_HEIGHT * 4
    width = panel_width * len(panels)

    body = []
    for i, (label, loci, length) in enumerate(panels):
        cx = i * panel_width + panel_width / 2
        cy = RING_RADIUS + PANEL_MARGIN
        body.append(_render_panel(cx, cy, loci, length, label))
    body.append(_render_legend(10, height - LEGEND_ROW_HEIGHT * 4 + 10))

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {width} {height}" width="100%">'
        f'{"".join(body)}</svg>'
    )


# ── entry point ──────────────────────────────────────────────────────

def build_panels(gff_path, annotation_summary_path, bin_metadata_path):
    cds_scores = load_cds_scores(annotation_summary_path)
    metadata = load_bin_metadata(bin_metadata_path)

    other_feats = parse_non_cds_features(gff_path)
    cds_feats = cds_features(gff_path, cds_scores)
    loci = cluster_loci(cds_feats, other_feats)

    length = genome_length(metadata, loci)
    if not length:
        return []

    panels = [('path1', loci, length)]
    segments = plastid_path2_segments(metadata)
    if segments is not None:
        panels.append(('path2', remap_to_path2(loci, segments), length))
    return panels


def run(gff_path, annotation_summary_path, bin_metadata_path, out_path):
    panels = build_panels(
        gff_path, annotation_summary_path, bin_metadata_path)
    if not panels:
        Path(out_path).write_text('')
        return
    Path(out_path).write_text(render_svg(panels))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff', required=True, type=Path)
    parser.add_argument('--annotation-summary', type=Path)
    parser.add_argument('--bin-metadata', type=Path)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()

    try:
        run(args.gff, args.annotation_summary, args.bin_metadata, args.out)
    except Exception as exc:  # noqa: BLE001 -- diagnostic must not fail
        print(
            f'organelle_map: {args.sample_id}: {exc}', file=sys.stderr)
        args.out.write_text('')

    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
