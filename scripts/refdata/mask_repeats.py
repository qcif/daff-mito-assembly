#!/usr/bin/env python3
"""
Mask low-complexity and microsatellite regions out of a recruit panel.

Tandem-repeat tails in some mitogenome references recruit nuclear reads
carrying the same repeat, inflating the recruited base count the
coverage gate sees. This replaces the repeat regions reported by
dustmasker with `N`. Length is preserved, so coordinates are unchanged.

Deliberately tool-free: the dustmasker call lives in the caller
(scripts/refdata/build_recruit.sh). Reuses the masking core of
mask_panel.py.

Usage:
    python3 mask_repeats.py \
        --fasta   animal_mt.raw.fa \
        --dust    animal_mt.dust.intervals \
        --out     animal_mt.fa \
        --report  animal_mt.masking.json \
        [--dust-level 30] [--dust-window 64] [--dust-linker 1] \
        [--max-masked-frac 0.60]
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mask_panel import (  # noqa: E402
    DEFAULT_MAX_MASKED_FRAC, mask_fasta, record_id)

# Recorded in the report only; the values passed to dustmasker live in
# build_recruit.sh.
DEFAULT_DUST_LEVEL = 30
DEFAULT_DUST_WINDOW = 64
DEFAULT_DUST_LINKER = 1

DUST_RANGE_SEP = " - "


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True, help="panel FASTA to mask")
    p.add_argument("--dust", required=True,
                   help="dustmasker -outfmt interval output for --fasta")
    p.add_argument("--out", required=True, help="masked FASTA output")
    p.add_argument("--report", required=True, help="JSON report output")
    p.add_argument("--dust-level", type=int, default=DEFAULT_DUST_LEVEL)
    p.add_argument("--dust-window", type=int, default=DEFAULT_DUST_WINDOW)
    p.add_argument("--dust-linker", type=int, default=DEFAULT_DUST_LINKER)
    p.add_argument("--max-masked-frac", type=float,
                   default=DEFAULT_MAX_MASKED_FRAC)
    return p.parse_args()


def load_dust_intervals(path: Path) -> dict:
    """
    Return {accession: [(start, end), ...]} from dustmasker interval output.

    The format is a `>header` line followed by `start - end` lines,
    0-based and inclusive of `end`; they are converted to half-open.
    """
    intervals = {}
    acc = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                acc = record_id(line)
                continue
            start, _, end = line.partition(DUST_RANGE_SEP)
            intervals.setdefault(acc, []).append((int(start), int(end) + 1))
    return intervals


def mask_repeats(
    fasta: Path,
    dust: Path,
    out_fasta: Path,
    max_masked_frac: float,
    dust_params: dict,
) -> dict:
    """Write a repeat-masked copy of `fasta`; return the masking report."""
    hits = load_dust_intervals(dust)
    report = mask_fasta(
        fasta, hits, out_fasta, max_masked_frac, label="repeat-masked")
    return {'masker': 'dustmasker', **dust_params, **report}


def main():
    args = parse_args()
    report = mask_repeats(
        fasta=Path(args.fasta),
        dust=Path(args.dust),
        out_fasta=Path(args.out),
        max_masked_frac=args.max_masked_frac,
        dust_params={
            'dust_level': args.dust_level,
            'dust_window': args.dust_window,
            'dust_linker': args.dust_linker,
        },
    )
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"[mask_repeats] {report['genomes_written']} records written, "
        f"{report['genomes_skipped']} skipped; "
        f"{report['bases_masked']} of {report['bases_written']} bases "
        f"masked ({report['masked_fraction']:.1%})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
