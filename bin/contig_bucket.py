"""Shared contig-classification -> bucket -> colour mapping (task 52).

Three views of a sample's binning decision -- the Assembly tab's
per-contig coverage chart (``bin/report/report.py``), its contig
table, and the assembly graph's node colouring (``annotate_graph_
svg.py``) -- all need to agree on what bucket a `bin_metadata.json`
classification falls into and what colour represents it. This module
is the single definition (CONSTITUTION rule 19); it is a plain
top-level ``bin/`` module rather than living under ``bin/report/``
because ``annotate_graph_svg.py`` runs in a container with no
``jinja2``/report dependencies and must not import anything that
pulls them in.

``MIXED_COLOUR`` is kept out of ``BUCKET_COLOURS`` deliberately: it is
a graph-only state (a GFA segment traversed by contigs from more than
one bucket), not a contig classification, so a contig-level consumer
like the coverage chart can never pick it up by accident.
"""

import json
from pathlib import Path
from typing import Optional

BUCKET_ORDER = ['target', 'secondary', 'off-target']

BUCKET_LABELS = {
    'target': 'Target',
    'secondary': 'Secondary',
    'off-target': 'Off-target',
}

BUCKET_COLOURS = {
    'target': '#2ca02c',
    'secondary': '#ff7f0e',
    'off-target': '#7f7f7f',
}

# Distinct from every bucket colour -- marks a graph segment traversed
# by contigs from more than one bucket (task 52 §3.4). Continues the
# tab10 palette the other three colours come from.
MIXED_COLOUR = '#9467bd'


def contig_bucket(classification: Optional[str]) -> str:
    """Maps a `bin_metadata.json` contig `classification` to one of
    the three display buckets. Anything other than the two explicit
    "the pipeline picked this" classifications -- including
    `off_target`, `sibling_organelle`, `None` and any future/unknown
    value -- reads as `off-target` (spec §3.7: those are all "not the
    target", and the graph/chart/table must agree on that)."""
    if classification == 'target_candidate':
        return 'target'
    if classification == 'secondary_target':
        return 'secondary'
    return 'off-target'


def bucket_legend() -> list:
    """`{key, name, color}` rows in display order -- feeds both the
    coverage chart (JS) and the Assembly tab's static colour key."""
    return [
        {'key': key, 'name': BUCKET_LABELS[key], 'color': BUCKET_COLOURS[key]}
        for key in BUCKET_ORDER
    ]


def load_bin_metadata(bin_metadata_path) -> dict:
    """Missing / empty / malformed `bin_metadata.json` -> `{}` --
    same degrade shape as `organelle_map.load_bin_metadata()`."""
    if not bin_metadata_path:
        return {}
    path = Path(bin_metadata_path)
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}
