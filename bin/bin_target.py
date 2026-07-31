#!/usr/bin/env python3
"""BIN_TARGET (C3) — spec §2 stage 10, §3.3.

Selects the dominant target contig from a METAFLYE assembly using the
intersection of:
  (a) coverage spike vs. per-contig median (from assembly_info.txt),
  (b) reference identity + aligned fraction to ${target}.mmi (minimap2),
  (c) ORF integrity under target-appropriate NCBI genetic-code table(s).

Emits target.fasta, secondaries.tsv, and bin_metadata.json describing
the decision. animal_mt: also runs the end-overlap circularity check.
"""

import argparse
import json
import shutil
import statistics
import sys
from pathlib import Path

import mappy
from Bio import SeqIO
from Bio.Seq import Seq


MIN_REF_IDENTITY = 80.0
MIN_ALIGNED_FRAC = 0.30
COVERAGE_SPIKE_FACTOR = 2.0
MIN_ORF_AA = 100
END_OVERLAP_WINDOW = 300
END_OVERLAP_MIN_IDENTITY = 0.90
END_OVERLAP_MIN_BASES = 200

SECONDARIES_COLS = [
    'contig_id', 'length_bp', 'coverage',
    'ref_identity_pct', 'orf_ok', 'classification',
]


def parse_assembly_info(path: Path) -> dict:
    """Parse Flye assembly_info.txt → {contig_id: {length, coverage}}."""
    result = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            result[parts[0]] = {
                'length': int(parts[1]),
                'coverage': float(parts[2]),
            }
    return result


def align_to_ref(seq: str, aligner: mappy.Aligner) -> tuple[float, float]:
    """Return (best_identity_pct, aligned_fraction) for seq vs loaded ref."""
    hits = list(aligner.map(seq))
    if not hits:
        return 0.0, 0.0
    best = max(hits, key=lambda h: h.mlen)
    identity = 100.0 * best.mlen / best.blen if best.blen else 0.0
    aligned_frac = best.blen / len(seq) if len(seq) else 0.0
    return identity, aligned_frac


def longest_orf_aa(seq: str, tables: list) -> tuple[int, int]:
    """
    Return (best_orf_aa_len, best_table) across all 6 frames and tables.
    Counts amino acids between stop codons (* in translation output).
    """
    best_len = 0
    best_table = tables[0]
    bio_seq = Seq(seq)
    for table in tables:
        for strand in (bio_seq, bio_seq.reverse_complement()):
            for frame in range(3):
                fragment = strand[frame:]
                trimmed = fragment[:len(fragment) - len(fragment) % 3]
                if not trimmed:
                    continue
                aa = str(trimmed.translate(table=table, to_stop=False))
                current = 0
                for aa_char in aa:
                    if aa_char == '*':
                        if current > best_len:
                            best_len = current
                            best_table = table
                        current = 0
                    else:
                        current += 1
                if current > best_len:
                    best_len = current
                    best_table = table
    return best_len, best_table


def classify_contigs(
    info: dict,
    sequences: dict,
    aligner: mappy.Aligner,
    tables: list,
) -> list:
    """Compute per-contig signals and return enriched row dicts."""
    coverages = [v['coverage'] for v in info.values()]
    median_cov = statistics.median(coverages) if coverages else 0.0

    rows = []
    for cid, meta in info.items():
        seq = sequences.get(cid, '')
        identity, aligned_frac = align_to_ref(seq, aligner)
        orf_len, best_table = longest_orf_aa(seq, tables)
        orf_ok = orf_len >= MIN_ORF_AA

        passes_identity = (
            identity >= MIN_REF_IDENTITY
            and aligned_frac >= MIN_ALIGNED_FRAC
        )
        passes_coverage = (
            median_cov > 0
            and meta['coverage'] >= COVERAGE_SPIKE_FACTOR * median_cov
        )

        if not passes_identity:
            classification = 'off_target'
        elif not passes_coverage:
            classification = 'secondary_target'
        elif not orf_ok:
            classification = 'low_confidence'
        else:
            classification = 'target_candidate'

        rows.append({
            'contig_id': cid,
            'length_bp': meta['length'],
            'coverage': meta['coverage'],
            'ref_identity_pct': round(identity, 2),
            'aligned_frac': round(aligned_frac, 4),
            'orf_ok': orf_ok,
            'orf_aa_len': orf_len,
            'selected_genetic_code': best_table,
            'classification': classification,
        })
    return rows


def select_primary(rows: list, assembly_target: str) -> tuple[list, list]:
    """
    Pure function: split rows into (primaries, secondaries).

    plant_mt: all target_candidates are primary.
    animal_mt / plant_pt: single best by (coverage, identity, orf_length).
    Non-selected candidates are reclassified as secondary_target in secondaries.
    """
    candidates = [r for r in rows if r['classification'] == 'target_candidate']
    non_candidates = [r for r in rows if r['classification'] != 'target_candidate']

    if not candidates:
        return [], rows

    if assembly_target == 'plant_mt':
        return candidates, non_candidates

    best = max(
        candidates,
        key=lambda r: (r['coverage'], r['ref_identity_pct'], r['orf_aa_len']),
    )
    runners_up = [
        {**r, 'classification': 'secondary_target'}
        for r in candidates
        if r['contig_id'] != best['contig_id']
    ]
    return [best], non_candidates + runners_up


def check_circularity(seq: str) -> bool:
    """
    End-overlap circularity check: self-align first vs last END_OVERLAP_WINDOW bp.
    Returns True if overlap detected (identity ≥ 90 %, aligned bases ≥ 200).
    """
    if len(seq) < 2 * END_OVERLAP_WINDOW:
        return False
    prefix = seq[:END_OVERLAP_WINDOW]
    suffix = seq[-END_OVERLAP_WINDOW:]
    aligner = mappy.Aligner(seq=suffix, preset='map-ont', best_n=1)
    hits = list(aligner.map(prefix))
    if not hits:
        return False
    best = max(hits, key=lambda h: h.mlen)
    identity = best.mlen / best.blen if best.blen else 0.0
    return identity >= END_OVERLAP_MIN_IDENTITY and best.blen >= END_OVERLAP_MIN_BASES


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assembly', type=Path, required=True)
    p.add_argument('--assembly-info', type=Path, required=True)
    p.add_argument('--gfa', type=Path, required=True)
    p.add_argument('--organelle-ref', type=Path, required=True)
    p.add_argument('--sample-id', required=True)
    p.add_argument('--assembly-target', required=True)
    p.add_argument('--genetic-codes', required=True,
                   help='Comma-separated NCBI genetic-code table IDs')
    p.add_argument('--out-target', type=Path, required=True)
    p.add_argument('--out-secondaries', type=Path, required=True)
    p.add_argument('--out-metadata', type=Path, required=True)
    args = p.parse_args()

    tables = [int(t) for t in args.genetic_codes.split(',')]

    info = parse_assembly_info(args.assembly_info)
    sequences = {
        r.id: str(r.seq)
        for r in SeqIO.parse(args.assembly, 'fasta')
    }

    aligner = mappy.Aligner(str(args.organelle_ref), preset='map-ont', best_n=5)
    if not aligner:
        print(
            f'WARNING: failed to load reference {args.organelle_ref}',
            file=sys.stderr,
        )

    rows = classify_contigs(info, sequences, aligner, tables)
    primaries, secondaries = select_primary(rows, args.assembly_target)

    with open(args.out_target, 'w') as fh:
        for row in primaries:
            cid = row['contig_id']
            seq = sequences.get(cid, '')
            fh.write(f'>{cid}\n{seq}\n')

    with open(args.out_secondaries, 'w') as fh:
        fh.write('\t'.join(SECONDARIES_COLS) + '\n')
        for row in secondaries:
            fh.write('\t'.join(str(row[c]) for c in SECONDARIES_COLS) + '\n')

    metadata = {
        'sample_id': args.sample_id,
        'assembly_target': args.assembly_target,
        'n_contigs_total': len(rows),
        'n_target_selected': len(primaries),
        'contigs_selected': [r['contig_id'] for r in primaries],
    }
    if primaries:
        best = primaries[0]
        metadata['selected_genetic_code'] = best['selected_genetic_code']
        metadata['ref_identity_pct'] = best['ref_identity_pct']
        metadata['coverage'] = best['coverage']
        metadata['orf_aa_len'] = best['orf_aa_len']

    if args.assembly_target == 'animal_mt' and primaries:
        primary_seq = sequences.get(primaries[0]['contig_id'], '')
        metadata['circular'] = check_circularity(primary_seq)

    if args.assembly_target == 'plant_pt':
        # Local import: keeps non-plant_pt code paths free of the C4
        # dependency (spec/plastid-canonicalisation.md §1, task 20 §2).
        from plastid_canonicalise import canonicalise_plastid

        iso_dir = args.out_target.resolve().parent / 'plastid_isoforms'
        result = canonicalise_plastid(args.gfa, outdir=iso_dir)

        metadata['plastid_canonicalisation'] = {
            **result._asdict(),
            'substitution_applied': result.branch == 'canonical',
        }

        if result.branch == 'canonical':
            shutil.copyfile(iso_dir / 'path1.fasta', args.out_target)
        elif iso_dir.exists():
            # C4 writes nothing off the canonical branch; a populated
            # directory here indicates something went wrong upstream.
            shutil.rmtree(iso_dir)
    else:
        metadata['plastid_canonicalisation'] = {
            'branch': 'not_applicable',
            'reason': f'assembly_target={args.assembly_target}',
        }

    args.out_metadata.write_text(json.dumps(metadata, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
