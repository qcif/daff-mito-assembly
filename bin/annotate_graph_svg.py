#!/usr/bin/env python3
"""Stage 9 -- C12: node identity for BandageNG's assembly graph SVG
(task 47).

BandageNG's ``image`` subcommand can write SVG, but the SVG it writes
carries no node identity at all: a Qt painter dump of anonymous
``<path>`` elements, no ``id`` attributes, no ``<text>``. There is also
no ``--names``/``--labels`` flag to add any. The one identity channel
BandageNG does expose is node colour: given a ``name,color`` CSV it
renders each GFA segment as exactly one path carrying that exact fill.

This module round-trips node identity through that channel via two
entry points, sharing one sentinel-allocation function so a rerun on
unchanged input is byte-identical (deterministic, derived from GFA
segment order -- never random):

- ``allocate`` -- GFA in, ``name,color`` sentinel CSV out. Runs before
  ``BandageNG image --color sentinels.csv``.
- ``annotate`` -- the rendered SVG + that CSV + ``assembly_info.txt``
  in; for each sentinel fill found, rewrites it to the neutral display
  colour and attaches ``data-node``, ``data-contigs`` and a ``<title>``
  child (a native SVG tooltip -- no JavaScript). Runs after BandageNG.

Both are their own Nextflow process (``ALLOCATE_GRAPH_SENTINELS`` /
``ANNOTATE_GRAPH_SVG``) either side of the unchanged ``BandageNG
image`` call, because the BandageNG biocontainer has no Python either
side of the render -- the same constraint that forced C9
(``select_genetic_code.py``) into its own step (CONSTITUTION rule 14).

**Never invents a per-node bin colour** -- a GFA segment (``edge_N``)
and a Flye contig (``contig_N``) are different namespaces joined
many-to-many through ``assembly_info.txt``'s ``graph_path`` column,
so a repeat edge shared between a target and an off-target contig has
no single correct bucket. Every traversing contig is listed in
``data-contigs`` instead, and the colouring stays neutral
(CONSTITUTION principle 7).

``annotate`` never raises past ``main()`` -- the graph is a diagnostic
and must not abort a sample (CONSTITUTION principle 8, rule 8). On any
failure it leaves the input SVG untouched rather than deleting it or
emitting nothing. A sentinel with no matching path, or a path with no
matching sentinel, are both warnings to stderr, not errors.
"""

import argparse
import re
import sys
from pathlib import Path
from xml.sax.saxutils import escape

# Neutral display colour every sentinel fill is rewritten to -- a mid
# grey that reads clearly against BandageNG's white background/black
# strokes without implying any bin/classification meaning.
DISPLAY_COLOUR = '#c8c8c8'

# Sentinel colours are allocated as #000001, #000002, ... (1-based, so
# 0 -- pure black, BandageNG's stroke colour -- is never assigned) and
# must never collide with the background (#ffffff) or link strokes
# (#000000); #ffffff is unreachable short of ~16.7M segments.
_SENTINEL_FLOOR = 1

# BandageNG's Qt SVG export puts each node's fill on the <g> element
# that wraps its <path> (and any siblings), not on the <path> itself
# (confirmed against the pinned biocontainer -- task 47 §2/§4, the
# brief's own before/after example is illustrative only, not literal
# output). Background/stroke elements never carry a `fill`, so this
# only ever matches node groups.
_FILL_RE = re.compile(r'<g\b[^>]*?\bfill="(#[0-9a-fA-F]{6})"')
_G_OPEN_TAG_RE = re.compile(r'<g\b[^>]*?\bfill="#[0-9a-fA-F]{6}"[^>]*>')


# ── allocate ─────────────────────────────────────────────────────────

def parse_gfa_segments(gfa_path) -> list:
    """Segment names in file order -- GFA `S` records, tab-separated,
    column 2. Order is the determinism anchor for sentinel allocation."""
    segments = []
    with open(gfa_path) as fh:
        for raw_line in fh:
            fields = raw_line.rstrip('\n').split('\t')
            if len(fields) < 2 or fields[0] != 'S':
                continue
            segments.append(fields[1])
    return segments


def allocate_sentinels(segments: list) -> list:
    """One unique, deterministic ``#RRGGBB`` per segment, in input
    order -- returns a list of ``(name, colour)`` pairs."""
    return [
        (name, f'#{i:06x}')
        for i, name in enumerate(segments, start=_SENTINEL_FLOOR)
    ]


def write_sentinel_csv(sentinels: list, out_path) -> None:
    """BandageNG rejects a headerless colour CSV (`didn't contain
    color`), so the header row is mandatory, not cosmetic."""
    lines = ['name,color']
    lines.extend(f'{name},{colour}' for name, colour in sentinels)
    Path(out_path).write_text('\n'.join(lines) + '\n')


def run_allocate(gfa_path, out_path) -> None:
    segments = parse_gfa_segments(gfa_path)
    write_sentinel_csv(allocate_sentinels(segments), out_path)


# ── contig traversal (assembly_info.txt graph_path) ────────────────

def parse_graph_path_column(assembly_info_path) -> dict:
    """Maps segment name -> sorted set of contig names that traverse
    it, from `assembly_info.txt` column 8 (`graph_path`).

    `graph_path` entries are signed edge numbers (`-2,1,2`); the sign
    marks traversal direction, not identity, so it is stripped (`-2`
    and `2` are the same segment `edge_2`). `*` means "no edge" and is
    skipped -- seen in real fixture data (`INT-ANIMAL-01`: `*,5,*`)."""
    traversal = {}
    if not assembly_info_path or not Path(assembly_info_path).is_file():
        return traversal
    with open(assembly_info_path) as fh:
        for raw_line in fh:
            if raw_line.startswith('#'):
                continue
            fields = raw_line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue
            contig_name = fields[0]
            for entry in fields[7].split(','):
                entry = entry.strip()
                if entry in ('', '*'):
                    continue
                try:
                    edge_num = int(entry.lstrip('-+'))
                except ValueError:
                    continue
                segment = f'edge_{edge_num}'
                traversal.setdefault(segment, set()).add(contig_name)
    return traversal


# ── annotate ─────────────────────────────────────────────────────────

_ATTR_ENTITIES = {'"': '&quot;'}

# BandageNG's raw SVG opens with an XML declaration (`<?xml version=
# "1.0" ...?>`) -- confirmed against the pinned biocontainer. That is
# fine as a standalone file but breaks when inlined mid-document (task
# 47 §4.1), so it -- and any DOCTYPE, which BandageNG does not emit
# but a future version might -- is stripped before publishing.
_XML_PROLOG_RE = re.compile(
    r'^\s*(?:<\?xml\b[^>]*\?>\s*)?(?:<!DOCTYPE\b[^>]*>\s*)?', re.I)


def strip_xml_prolog(svg_text: str) -> str:
    return _XML_PROLOG_RE.sub('', svg_text, count=1)


def _recoloured_open_tag(
    open_tag: str, colour: str, name: str, contigs: set,
) -> str:
    contig_names = ','.join(sorted(contigs))
    contig_attr = escape(contig_names, _ATTR_ENTITIES)
    node_attr = escape(name, _ATTR_ENTITIES)
    title = escape(name if not contigs else f'{name} — {contig_names}')
    recoloured = open_tag.replace(
        f'fill="{colour}"', f'fill="{DISPLAY_COLOUR}"', 1)
    # Inject the identity attributes into the <g ...> opening tag, and
    # a <title> child immediately after it -- SVG renders <title> as a
    # native browser tooltip on its parent element with no JS involved.
    tagged = recoloured[:-1] + (
        f' data-node="{node_attr}" data-contigs="{contig_attr}">')
    return f'{tagged}<title>{title}</title>'


def annotate_svg(
    svg_text: str, sentinels: list, contig_traversal: dict,
    warn=None,
) -> str:
    """Rewrites each sentinel-coloured ``<g fill="#rrggbb" ...>`` node
    wrapper in place: the fill is replaced with the neutral display
    colour, ``data-node``/``data-contigs`` attributes are added, and a
    ``<title>`` child is inserted as the wrapper's first child.
    Unmatched sentinels/fills are reported via ``warn`` (defaults to a
    no-op) rather than raised."""
    if warn is None:
        def warn(_msg):
            return None

    colour_to_name = {colour: name for name, colour in sentinels}
    fills_in_svg = set(_FILL_RE.findall(svg_text))

    for colour, name in ((c, n) for n, c in sentinels):
        if colour not in fills_in_svg:
            warn(f'sentinel {colour} ({name}) matches no path in the SVG')

    def _replace_group(match: re.Match) -> str:
        open_tag = match.group(0)
        fill_match = re.search(r'fill="(#[0-9a-fA-F]{6})"', open_tag)
        colour = fill_match.group(1)
        name = colour_to_name.get(colour)
        if name is None:
            if colour not in ('#ffffff', '#000000'):
                warn(f'path fill {colour} matches no sentinel')
            return open_tag
        return _recoloured_open_tag(
            open_tag, colour, name, contig_traversal.get(name, set()))

    return _G_OPEN_TAG_RE.sub(_replace_group, svg_text)


def run_annotate(
    svg_path, sentinels_csv_path, assembly_info_path, out_path,
    sample_id: str,
) -> None:
    svg_text = strip_xml_prolog(Path(svg_path).read_text())
    sentinels = read_sentinel_csv(sentinels_csv_path)
    contig_traversal = parse_graph_path_column(assembly_info_path)

    def warn(msg):
        print(f'annotate_graph_svg: {sample_id}: {msg}', file=sys.stderr)

    annotated = annotate_svg(svg_text, sentinels, contig_traversal, warn)
    Path(out_path).write_text(annotated)


def read_sentinel_csv(path) -> list:
    sentinels = []
    with open(path) as fh:
        next(fh, None)  # header
        for raw_line in fh:
            fields = raw_line.rstrip('\n').split(',')
            if len(fields) != 2:
                continue
            sentinels.append((fields[0], fields[1]))
    return sentinels


# ── entry point ──────────────────────────────────────────────────────

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)

    p_allocate = sub.add_parser(
        'allocate', help='GFA in, name,color sentinel CSV out')
    p_allocate.add_argument('--gfa', required=True, type=Path)
    p_allocate.add_argument('--out', required=True, type=Path)

    p_annotate = sub.add_parser(
        'annotate',
        help='SVG + sentinel CSV + assembly_info.txt in, annotated SVG out')
    p_annotate.add_argument('--svg', required=True, type=Path)
    p_annotate.add_argument('--sentinels', required=True, type=Path)
    p_annotate.add_argument('--assembly-info', type=Path, default=None)
    p_annotate.add_argument('--sample-id', required=True)
    p_annotate.add_argument('--out', required=True, type=Path)

    args = parser.parse_args()

    if args.command == 'allocate':
        run_allocate(args.gfa, args.out)
        return 0

    # annotate: always exits 0 -- a rendering defect must not fail the
    # sample (CONSTITUTION principle 8). On any failure, leave the
    # un-annotated input SVG in place rather than deleting or emitting
    # nothing (rule 18).
    try:
        run_annotate(
            args.svg, args.sentinels, args.assembly_info, args.out,
            args.sample_id,
        )
    except Exception as exc:  # noqa: BLE001 -- diagnostic must not fail
        print(
            f'annotate_graph_svg: {args.sample_id}: {exc}', file=sys.stderr)
        if Path(args.svg).is_file():
            args.out.write_text(Path(args.svg).read_text())
        else:
            args.out.write_text('')

    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
