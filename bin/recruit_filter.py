#!/usr/bin/env python3
"""RECRUIT read selection — spec §2 stage 5, §9 item 1.

Reads a PAF of cleaned reads aligned against the declared organelle
panel and emits the IDs of reads to recruit, one per line.

Selection is on **merged aligned extent**, never on mapping quality.
RECRUIT previously filtered with `samtools view -F 4 -q 1`, which
discards MAPQ 0. That worked while the panel held a single reference
genome — MAPQ 0 then meant "this read came from a repeat, I cannot say
where" — but it inverts on the broad, redundant panels task 28
introduced. minimap2 derives MAPQ from the gap between the best and
second-best chaining score, so on 652 near-identical mitogenomes the
runner-up is *the same locus in another species*. A read matching
conserved organellar sequence present in hundreds of panel genomes
scores MAPQ 0; a plastid read matching residual unmasked NUPT in one
genome scores MAPQ 60. Measured on `INT-PLANT-01-mt`, the MAPQ-0 pool
was 54 % mitochondrial and the MAPQ-60 pool only 8.7 % — the filter
selected against the target organelle (task 28 §9.5).

Merged aligned extent asks the question recruitment actually cares
about — *how much of this read looks like the target organelle* — and
its answer does not change when the panel gains more genomes.

Both thresholds are per-target and come from
`params.recruit_thresholds`, not from constants here
(CONSTITUTION rule 18). Defaults of 0 admit every read with any
alignment, which is the measured-correct setting today; see §9 item 1
for the benchmark that should set them.
"""

import argparse
import sys
from pathlib import Path

# Sibling module in bin/, staged as a directory onto the container PATH.
from intervals import merged_length


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--paf", required=True,
                   help="reads aligned against the declared panel")
    p.add_argument("--out", default="-",
                   help="read-ID list output ('-' for stdout)")
    p.add_argument("--min-aligned-frac", type=float, default=0.0,
                   help="minimum merged aligned bases / read length")
    p.add_argument("--min-aligned-bp", type=int, default=0,
                   help="minimum merged aligned bases")
    p.add_argument("--stats", default="",
                   help="optional JSON path for selection counts")
    return p.parse_args()


def read_alignments(paf: Path) -> dict:
    """
    Return {read_id: (read_length, merged_aligned_bp)} from a PAF.

    PAF columns 1/2/3/4 are query name, query length, query start and
    query end (0-based, half-open). Query intervals are merged across
    every alignment block for the read, so a read matching the panel in
    several places is counted once — the same metric C2 and C3 use
    (spec §3.7.2).
    """
    lengths, blocks = {}, {}
    with open(paf) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.split('\t')
            read_id = fields[0]
            lengths[read_id] = int(fields[1])
            blocks.setdefault(read_id, []).append(
                (int(fields[2]), int(fields[3])))
    return {
        read_id: (lengths[read_id], merged_length(spans))
        for read_id, spans in blocks.items()
    }


def select(alignments: dict, min_frac: float, min_bp: int) -> list:
    """Read IDs clearing both thresholds, sorted for stable output."""
    keep = []
    for read_id, (length, aligned) in alignments.items():
        if aligned < min_bp:
            continue
        if length and aligned / length < min_frac:
            continue
        keep.append(read_id)
    return sorted(keep)


def main():
    args = parse_args()
    alignments = read_alignments(Path(args.paf))
    keep = select(alignments, args.min_aligned_frac, args.min_aligned_bp)

    handle = sys.stdout if args.out == "-" else open(args.out, 'w')
    try:
        handle.write("".join(f"{read_id}\n" for read_id in keep))
    finally:
        if handle is not sys.stdout:
            handle.close()

    print(
        f"[recruit_filter] {len(keep)} of {len(alignments)} aligned reads "
        f"recruited (min_aligned_frac={args.min_aligned_frac}, "
        f"min_aligned_bp={args.min_aligned_bp})",
        file=sys.stderr,
    )

    if args.stats:
        import json
        Path(args.stats).write_text(json.dumps({
            'reads_aligned': len(alignments),
            'reads_recruited': len(keep),
            'min_aligned_frac': args.min_aligned_frac,
            'min_aligned_bp': args.min_aligned_bp,
        }, indent=2) + "\n")


if __name__ == "__main__":  # pragma: no cover
    main()
