#!/usr/bin/env python3
"""Plastid quadripartite canonicalisation (C4).

See spec/plastid-canonicalisation.md.

Reads a Flye ``assembly_graph.gfa`` and classifies the plant_pt plastid
graph into one of three branches based on its S-line edge count:

  * 1 edge  -> "resolved_circle" (already a single closed plastome).
  * 3 edges -> "canonical" if the length/depth degenerate checks pass;
               writes path1.fasta (LSC+IR+SSC+rc(IR)) and path2.fasta
               (LSC+IR+rc(SSC)+rc(IR)), the two biologically valid SSC
               orientations.
  * anything else -> "non_canonical", diagnostic only.

Invoked in-process by bin/bin_target.py (C3) on the plant_pt branch of
stage 10 (BIN_TARGET); see spec/plastid-canonicalisation.md for the
full algorithm specification, which is the sole permitted reference
for this implementation.
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import NamedTuple, Optional

from Bio.Seq import Seq


# spec §3.3/§10: case-insensitive on tag name, type-agnostic on the GFA
# type character (real Flye output uses dp:i:, not dp:f:); no scientific
# notation.
DEPTH_TAG_RE = re.compile(r'[Dd][Pp]:[a-zA-Z]:([0-9.]+)')


class Edge(NamedTuple):
    name: str
    sequence: str
    depth: float
    length: int


class Result(NamedTuple):
    branch: str
    edge_count: int
    lsc_edge: Optional[str]
    ir_edge: Optional[str]
    ssc_edge: Optional[str]
    path1_len: Optional[int]
    path2_len: Optional[int]
    non_canonical_reason: Optional[str]


def parse_gfa(gfa_path: Path) -> list:
    """Parse S lines from a GFA file into a list of Edge tuples.

    Ignores L/H/P/W/etc. lines entirely (spec §3.1). A missing depth
    tag is treated as depth 0.0 (spec §3.3); a malformed depth value
    is also treated as 0.0 (spec §10 implementation choice).
    """
    edges = []
    seen_names = set()
    with open(gfa_path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip('\r\n')
            fields = line.split('\t')
            if not fields or fields[0] != 'S':
                continue

            name = fields[1]
            sequence = fields[2]
            if sequence == '*':
                raise ValueError(
                    f"edge {name!r}: '*' (absent sequence) is not supported"
                )
            if name in seen_names:
                raise ValueError(f"duplicate edge name in GFA: {name!r}")
            seen_names.add(name)

            depth = 0.0
            for tag in fields[3:]:
                match = DEPTH_TAG_RE.search(tag)
                if match:
                    try:
                        depth = float(match.group(1))
                    except ValueError:
                        depth = 0.0
                    break

            edges.append(
                Edge(name=name, sequence=sequence, depth=depth,
                     length=len(sequence))
            )
    return edges


def _select_lsc(edges: list) -> Edge:
    """Longest edge is LSC.

    Tie-break (spec §6.1, implementation choice): length desc, then
    depth desc, then edge name ascending (lexicographically first).
    """
    max_length = max(e.length for e in edges)
    candidates = [e for e in edges if e.length == max_length]
    if len(candidates) == 1:
        return candidates[0]
    max_depth = max(e.depth for e in candidates)
    candidates = [e for e in candidates if e.depth == max_depth]
    if len(candidates) == 1:
        return candidates[0]
    return min(candidates, key=lambda e: e.name)


def _select_ir(edges: list) -> Edge:
    """Deepest edge is IR.

    Tie-break (spec §6.1, implementation choice): depth desc, then
    length desc, then edge name descending (lexicographically last).
    """
    max_depth = max(e.depth for e in edges)
    candidates = [e for e in edges if e.depth == max_depth]
    if len(candidates) == 1:
        return candidates[0]
    max_length = max(e.length for e in candidates)
    candidates = [e for e in candidates if e.length == max_length]
    if len(candidates) == 1:
        return candidates[0]
    return max(candidates, key=lambda e: e.name)


def _non_canonical(edge_count: int, reason: str) -> Result:
    return Result(
        branch='non_canonical',
        edge_count=edge_count,
        lsc_edge=None,
        ir_edge=None,
        ssc_edge=None,
        path1_len=None,
        path2_len=None,
        non_canonical_reason=reason,
    )


def _classify_three_edge(edges: list, outdir: Path) -> Result:
    """Apply the 3-edge degenerate checks (spec §5.2) then canonicalise."""
    # Check 1: depth collision — deepest edge must be strictly deeper
    # than at least one other edge.
    if len({e.depth for e in edges}) == 1:
        return _non_canonical(
            3, 'depth_tie: all edges have equal depth'
        )

    lsc = _select_lsc(edges)
    ir = _select_ir(edges)

    # Check 2: LSC and IR must not be the same edge.
    if lsc.name == ir.name:
        return _non_canonical(
            3, 'lsc_ir_collision: longest and deepest are the same edge'
        )

    # Check 3: zero-length edge is invalid input.
    if any(e.length == 0 for e in edges):
        return _non_canonical(3, 'zero_length_edge')

    ssc = next(e for e in edges if e.name not in (lsc.name, ir.name))

    ir_rc = str(Seq(ir.sequence).reverse_complement())
    ssc_rc = str(Seq(ssc.sequence).reverse_complement())
    path1_seq = lsc.sequence + ir.sequence + ssc.sequence + ir_rc
    path2_seq = lsc.sequence + ir.sequence + ssc_rc + ir_rc

    expected_len = lsc.length + 2 * ir.length + ssc.length
    if len(path1_seq) != expected_len or len(path2_seq) != expected_len:
        raise ValueError(
            'path length invariant failed: '
            f'path1={len(path1_seq)} path2={len(path2_seq)} '
            f'expected={expected_len}'
        )

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'path1.fasta').write_text(f'>path1\n{path1_seq}\n')
    (outdir / 'path2.fasta').write_text(f'>path2\n{path2_seq}\n')

    return Result(
        branch='canonical',
        edge_count=3,
        lsc_edge=lsc.name,
        ir_edge=ir.name,
        ssc_edge=ssc.name,
        path1_len=len(path1_seq),
        path2_len=len(path2_seq),
        non_canonical_reason=None,
    )


def canonicalise_plastid(gfa_path, outdir='.') -> Result:
    """Classify and (if canonical) canonicalise a plastid assembly graph.

    See spec/plastid-canonicalisation.md §8 for the full contract.
    """
    gfa_path = Path(gfa_path)
    if not gfa_path.exists():
        raise FileNotFoundError(str(gfa_path))
    outdir = Path(outdir)

    edges = parse_gfa(gfa_path)
    edge_count = len(edges)

    if edge_count == 1:
        return Result(
            branch='resolved_circle',
            edge_count=1,
            lsc_edge=None,
            ir_edge=None,
            ssc_edge=None,
            path1_len=None,
            path2_len=None,
            non_canonical_reason=None,
        )

    if edge_count != 3:
        return _non_canonical(edge_count, f'edge_count:{edge_count}')

    return _classify_three_edge(edges, outdir)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('gfa', type=Path, help='Path to assembly_graph.gfa')
    p.add_argument('--outdir', type=Path, default=Path('.'),
                   help='Directory for output FASTA files')
    p.add_argument('--json-out', type=Path, default=None,
                   help='Write the Result metadata as JSON to this path')
    args = p.parse_args()

    result = canonicalise_plastid(args.gfa, outdir=args.outdir)

    summary = f'branch={result.branch}'
    if result.branch == 'canonical':
        summary += (
            f' lsc={result.lsc_edge} ir={result.ir_edge}'
            f' ssc={result.ssc_edge}'
        )
    elif result.branch == 'non_canonical':
        summary += (
            f' edge_count={result.edge_count}'
            f' reason={result.non_canonical_reason}'
        )

    print(summary)

    if args.json_out:
        args.json_out.write_text(json.dumps(result._asdict(), indent=2))

    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
