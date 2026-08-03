#!/usr/bin/env python3
"""
Mask plastid-derived regions out of a plant mitochondrion panel.

Plant mitogenomes carry NUPTs — chunks of chloroplast DNA copied into
the mitochondrial genome over evolutionary time. Across RefSeq
Viridiplantae these account for a quarter of the total sequence, so a
broad, unmodified plant-mitochondrion recruitment panel is in part a
chloroplast panel: plastid contigs then score against the mitochondrial
reference as well as they do against the chloroplast one, and the
sibling-organelle discrimination in spec §3.7.1 collapses (task 28 §1).

This script keeps the breadth and removes the confusion. Given a PAF of
mitogenomes aligned against the plastid panel, it merges the aligned
query intervals per mitogenome (spec §3.7.2, the same routine C2 and C3
use — see bin/intervals.py) and replaces them with `N`.

Masking substitutes, never deletes: every output record has exactly the
length of its input, so coordinates are preserved.

Deliberately tool-free. The minimap2 call lives in the caller
(scripts/refdata/build_recruit.sh) so this stays pure-Python and its
tests need no boundary mocking (task 27 §6).

Usage:
    python3 mask_panel.py \
        --fasta   refseq_mt_viridiplantae.fa \
        --paf     mt_vs_pt.paf \
        --out     plant_mt.fa \
        --report  plant_mt.masking.json \
        [--min-align-len 200] \
        [--max-masked-frac 0.60]
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(
    0, str(Path(__file__).resolve().parents[2] / "bin"))

from intervals import merge_intervals, merged_length  # noqa: E402

# Output FASTA line width. Cosmetic — the panel is consumed by
# minimap2, which is agnostic to wrapping.
FASTA_LINE_WIDTH = 60

# Below this many aligned query bases an interval is not masked: short
# hits are conserved-gene fragments shared by both organelles rather
# than genuine NUPT insertions (task 28 §3.2).
DEFAULT_MIN_ALIGN_LEN = 200

# A mitogenome measuring above this fraction plastid is more likely a
# misannotated RefSeq record than a real genome. Skip it rather than
# emit a mostly-`N` reference. A data-quality guard, not a tuning knob.
DEFAULT_MAX_MASKED_FRAC = 0.60

# How many worst-masked genomes to name in the report.
TOP_OFFENDERS = 20

MASK_CHAR = 'N'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True,
                   help="mitogenome FASTA to mask")
    p.add_argument("--paf", required=True,
                   help="PAF of --fasta aligned against the plastid panel")
    p.add_argument("--out", required=True,
                   help="masked FASTA output path")
    p.add_argument("--report", required=True,
                   help="JSON masking report output path")
    p.add_argument("--min-align-len", type=int,
                   default=DEFAULT_MIN_ALIGN_LEN)
    p.add_argument("--max-masked-frac", type=float,
                   default=DEFAULT_MAX_MASKED_FRAC)
    return p.parse_args()


def iter_fasta_records(path):
    """Yield (header_line, sequence) with the sequence unwrapped."""
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def record_id(header: str) -> str:
    """First whitespace-delimited token of a FASTA header, sans '>'."""
    return header[1:].split()[0]


def load_mask_intervals(paf: Path, min_align_len: int) -> dict:
    """
    Return {query_id: [(q_start, q_end), ...]} from a PAF.

    PAF columns 1/3/4 are query name, query start and query end
    (0-based, half-open). Intervals spanning fewer than
    `min_align_len` query bases are dropped before merging.
    """
    intervals = {}
    with open(paf) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.split('\t')
            q_start, q_end = int(fields[2]), int(fields[3])
            if q_end - q_start < min_align_len:
                continue
            intervals.setdefault(fields[0], []).append((q_start, q_end))
    return intervals


def mask_sequence(seq: str, intervals: list) -> str:
    """
    Replace every merged interval of `seq` with `N`.

    Intervals are clamped to the sequence; the returned string always
    has len(seq) — masking substitutes, never deletes.
    """
    if not intervals:
        return seq
    chars = list(seq)
    for start, end in merge_intervals(intervals):
        start, end = max(0, start), min(len(seq), end)
        for i in range(start, end):
            chars[i] = MASK_CHAR
    return "".join(chars)


def clamped_masked_length(seq_len: int, intervals: list) -> int:
    """Merged length of `intervals` after clamping to [0, seq_len)."""
    clamped = [
        (max(0, s), min(seq_len, e))
        for s, e in intervals
        if max(0, s) < min(seq_len, e)
    ]
    return merged_length(clamped)


def wrap(seq: str, width: int = FASTA_LINE_WIDTH) -> str:
    """Wrap a sequence to `width` columns."""
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def mask_panel(
    fasta: Path,
    paf: Path,
    out_fasta: Path,
    min_align_len: int,
    max_masked_frac: float,
) -> dict:
    """
    Write a plastid-masked copy of `fasta`; return the masking report.

    Genomes measuring above `max_masked_frac` plastid are omitted from
    the output entirely and recorded under `skipped_genomes`.
    """
    hits = load_mask_intervals(paf, min_align_len)

    per_genome = []
    n_written = n_skipped = 0
    total_bases = total_masked = 0

    with open(out_fasta, 'w') as fh:
        for header, seq in iter_fasta_records(fasta):
            acc = record_id(header)
            intervals = hits.get(acc, [])
            masked_bp = clamped_masked_length(len(seq), intervals)
            frac = masked_bp / len(seq) if seq else 0.0
            per_genome.append({
                'accession': acc,
                'length_bp': len(seq),
                'masked_bp': masked_bp,
                'masked_fraction': round(frac, 6),
                'skipped': frac > max_masked_frac,
            })

            if frac > max_masked_frac:
                print(
                    f"[mask_panel] WARNING: skipping {acc} — "
                    f"{frac:.1%} plastid-derived, above the "
                    f"{max_masked_frac:.0%} guard",
                    file=sys.stderr,
                )
                n_skipped += 1
                continue

            fh.write(f"{header}\n{wrap(mask_sequence(seq, intervals))}\n")
            n_written += 1
            total_bases += len(seq)
            total_masked += masked_bp

    offenders = sorted(
        per_genome, key=lambda g: g['masked_fraction'], reverse=True)

    return {
        'min_align_len': min_align_len,
        'max_masked_frac': max_masked_frac,
        'genomes_in': len(per_genome),
        'genomes_written': n_written,
        'genomes_skipped': n_skipped,
        'genomes_with_plastid_hits': sum(
            1 for g in per_genome if g['masked_bp'] > 0),
        'bases_written': total_bases,
        'bases_masked': total_masked,
        'masked_fraction': round(
            total_masked / total_bases, 6) if total_bases else 0.0,
        'skipped_genomes': [g for g in per_genome if g['skipped']],
        'top_masked_genomes': offenders[:TOP_OFFENDERS],
    }


def main():
    args = parse_args()
    report = mask_panel(
        fasta=Path(args.fasta),
        paf=Path(args.paf),
        out_fasta=Path(args.out),
        min_align_len=args.min_align_len,
        max_masked_frac=args.max_masked_frac,
    )
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"[mask_panel] {report['genomes_written']} genomes written, "
        f"{report['genomes_skipped']} skipped; "
        f"{report['bases_masked']} of {report['bases_written']} bases "
        f"masked ({report['masked_fraction']:.1%})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
