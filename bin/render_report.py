#!/usr/bin/env python3
"""Standalone CLI for `report/report.py`'s `render()` (task 43a §5.1).

`collate.py` calls `render()` in-process (wrapped in its own
try/except — §5.2), so COLLATE never invokes this file. It exists as
the renderer's documented entry point for manual regeneration and for
task 45_run_report.md, which reuses the same Jinja machinery rather
than forking it.
"""

import argparse
import json
import sys
from pathlib import Path

from report.report import render


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata-json", type=Path, required=True)
    p.add_argument("--template-dir", type=Path, required=True)
    p.add_argument("--static-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--params-json", type=Path, default=None)
    p.add_argument("--nanoplot-raw", type=Path, default=None)
    p.add_argument("--nanoplot-clean", type=Path, default=None)
    p.add_argument("--organelle-map-svg", type=Path, default=None)
    p.add_argument("--graph-png", type=Path, default=None)
    p.add_argument("--annotation-gff", type=Path, default=None)
    args = p.parse_args()

    metadata = json.loads(args.metadata_json.read_text())
    params = (
        json.loads(args.params_json.read_text())
        if args.params_json and args.params_json.is_file() else {}
    )

    render(
        metadata=metadata,
        template_dir=args.template_dir,
        static_dir=args.static_dir,
        out_path=args.out,
        params=params,
        nanoplot_raw_dir=args.nanoplot_raw,
        nanoplot_clean_dir=args.nanoplot_clean,
        organelle_map_svg=args.organelle_map_svg,
        graph_png=args.graph_png,
        annotation_gff=args.annotation_gff,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
