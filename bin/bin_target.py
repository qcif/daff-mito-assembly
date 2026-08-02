#!/usr/bin/env python3
"""BIN_TARGET (C3) — spec §2 stage 10, §3.3, §3.7.

Selects the contig(s) constituting the declared organelle from a
METAFLYE assembly. A contig is a target candidate where all of
(spec §3.7.6):

  (a) merged aligned fraction against the declared panel
      >= --min-aligned-frac,
  (b) merged identity against the declared panel >= --min-identity,
  (c) merged aligned fraction against the declared panel exceeds that
      against every sibling organelle panel (spec §3.7.1).

Coverage is a ranking signal, not an admission gate (spec §3.7.3).
Longest-ORF length is recorded as a diagnostic but does not gate
selection (spec §3.7.5).

Emits target.fasta, secondaries.tsv, and bin_metadata.json describing
the decision, including the aligned fraction against every panel
scored, the winning panel, and how circularity was determined.
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

import mappy
from Bio import SeqIO
from Bio.Seq import Seq


# Sibling organelle panels scored alongside the declared panel, per
# spec §3.7.1. Plant samples carry both organelles; RECRUIT's positive
# selection does not separate them.
SIBLING_PANELS = {
    'plant_mt': ['plant_pt'],
    'plant_pt': ['plant_mt'],
    'animal_mt': [],
}

MIN_ORF_AA = 100
END_OVERLAP_WINDOW = 300
END_OVERLAP_MIN_IDENTITY = 0.90
END_OVERLAP_MIN_BASES = 200

FLYE_CIRC_YES = 'Y'

BASE_SECONDARIES_COLS = [
    'contig_id', 'length_bp', 'coverage',
    'ref_identity_pct', 'aligned_frac', 'winning_panel',
    'circular', 'circular_method', 'orf_aa_len',
    'low_coverage_candidate', 'classification',
]


def parse_assembly_info(path: Path) -> dict:
    """
    Parse Flye assembly_info.txt.

    Columns: seq_name, length, cov., circ., repeat, mult., ...
    Returns {contig_id: {length, coverage, flye_circular}}.
    """
    result = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            circ = parts[3].strip() if len(parts) > 3 else ''
            result[parts[0]] = {
                'length': int(parts[1]),
                'coverage': float(parts[2]),
                'flye_circular': circ.upper() == FLYE_CIRC_YES,
            }
    return result


def merge_intervals(intervals: list) -> list:
    """
    Merge overlapping and adjacent [start, end) intervals.

    Returns a sorted, non-overlapping list. Adjacent intervals
    (end == next start) are merged; empty input returns [].
    """
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def align_to_ref(seq: str, aligner: mappy.Aligner) -> tuple[float, float]:
    """
    Return (identity_pct, merged_aligned_fraction) for seq vs a panel.

    Aligned fraction is the union of query intervals over all alignment
    blocks divided by contig length; identity is aggregated as
    sum(matching bases) / sum(block lengths) over the same hits
    (spec §3.7.2).
    """
    if not seq:
        return 0.0, 0.0
    hits = list(aligner.map(seq))
    if not hits:
        return 0.0, 0.0
    merged = merge_intervals([(h.q_st, h.q_en) for h in hits])
    aligned_bp = sum(end - start for start, end in merged)
    total_blen = sum(h.blen for h in hits)
    total_mlen = sum(h.mlen for h in hits)
    identity = 100.0 * total_mlen / total_blen if total_blen else 0.0
    return identity, aligned_bp / len(seq)


def longest_orf_aa(seq: str, tables: list) -> tuple[int, int]:
    """
    Return (best_orf_aa_len, best_table) across all 6 frames and tables.
    Counts amino acids between stop codons (* in translation output).

    Reported provenance only — not a selection criterion (spec §3.7.5).
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


def check_circularity(seq: str) -> bool:
    """
    End-overlap circularity check: self-align the first vs the last
    END_OVERLAP_WINDOW bp. Returns True if an overlap is detected
    (identity >= 90 %, aligned bases >= 200).

    Fallback only — Flye trims the terminal overlap when it
    circularises, so this can fire only on contigs Flye marked 'N'
    (spec §3.7.4).
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
    return (
        identity >= END_OVERLAP_MIN_IDENTITY
        and best.blen >= END_OVERLAP_MIN_BASES
    )


def resolve_circularity(flye_circular: bool, seq: str) -> tuple[bool, str]:
    """Return (circular, method) per spec §3.7.4."""
    if flye_circular:
        return True, 'flye_circ'
    if check_circularity(seq):
        return True, 'end_overlap'
    return False, 'none'


def load_panel(mmi: Path):
    """
    Load one minimap2 index, or return None if it is unusable.

    mappy.Aligner is truthy even for a missing or corrupt index, so an
    empty `seq_names` is the reliable load-failure signal.
    """
    if not mmi.exists():
        return None
    aligner = mappy.Aligner(str(mmi), preset='map-ont', best_n=5)
    return aligner if aligner.seq_names else None


def load_panels(
    ref_dir: Path,
    assembly_target: str,
) -> tuple[dict, list]:
    """
    Load the declared panel plus its sibling panel(s) from the bundle.

    Returns (aligners, sibling_panels_scored). A missing sibling index
    warns to stderr and is skipped rather than failing the sample
    (CONSTITUTION principle 8). A missing declared index is fatal.
    """
    declared_mmi = ref_dir / f'{assembly_target}.mmi'
    declared = load_panel(declared_mmi)
    if declared is None:
        raise SystemExit(
            f'ERROR: failed to load declared panel {declared_mmi}')

    aligners = {assembly_target: declared}
    scored = []
    for sibling in SIBLING_PANELS.get(assembly_target, []):
        mmi = ref_dir / f'{sibling}.mmi'
        aligner = load_panel(mmi)
        if aligner is None:
            print(
                f'WARNING: sibling panel {mmi} is missing or unreadable'
                ' — sibling-organelle discrimination is disabled for'
                f' {assembly_target}',
                file=sys.stderr,
            )
            continue
        aligners[sibling] = aligner
        scored.append(sibling)
    return aligners, scored


def classify_contigs(
    info: dict,
    sequences: dict,
    aligners: dict,
    assembly_target: str,
    tables: list,
    thresholds: dict,
) -> list:
    """
    Compute per-contig signals and return enriched row dicts.

    `aligners` maps panel name → loaded mappy.Aligner and must contain
    the declared panel; any further entries are sibling panels.
    """
    siblings = [p for p in aligners if p != assembly_target]
    min_identity = thresholds['min_identity']
    min_aligned_frac = thresholds['min_aligned_frac']

    rows = []
    for cid, meta in info.items():
        seq = sequences.get(cid, '')
        panel_frac = {}
        panel_identity = {}
        for panel, aligner in aligners.items():
            identity, frac = align_to_ref(seq, aligner)
            panel_frac[panel] = round(frac, 4)
            panel_identity[panel] = round(identity, 2)

        declared_frac = panel_frac[assembly_target]
        declared_identity = panel_identity[assembly_target]
        best_sibling_frac = max(
            (panel_frac[p] for p in siblings), default=0.0)
        winning_panel = max(panel_frac, key=lambda p: panel_frac[p])
        if panel_frac[winning_panel] <= 0.0:
            winning_panel = 'none'

        orf_len, best_table = longest_orf_aa(seq, tables)
        circular, circular_method = resolve_circularity(
            meta['flye_circular'], seq)

        if best_sibling_frac > 0.0 and declared_frac <= best_sibling_frac:
            classification = 'sibling_organelle'
        elif (
            declared_frac < min_aligned_frac
            or declared_identity < min_identity
        ):
            classification = 'off_target'
        else:
            classification = 'target_candidate'

        rows.append({
            'contig_id': cid,
            'length_bp': meta['length'],
            'coverage': meta['coverage'],
            'ref_identity_pct': declared_identity,
            'aligned_frac': declared_frac,
            'panel_aligned_frac': panel_frac,
            'panel_identity_pct': panel_identity,
            'winning_panel': winning_panel,
            'circular': circular,
            'circular_method': circular_method,
            'orf_ok': orf_len >= MIN_ORF_AA,
            'orf_aa_len': orf_len,
            'selected_genetic_code': best_table,
            'low_coverage_candidate': False,
            'classification': classification,
        })
    return rows


def select_primary(
    rows: list,
    assembly_target: str,
    thresholds: dict,
) -> tuple[list, list]:
    """
    Pure function: split rows into (primaries, secondaries).

    emit == 'all' (plant_mt): every candidate up to max_contigs, ranked
    by coverage; candidates whose coverage falls below
    low_coverage_fraction × the top candidate's are emitted but flagged
    `low_coverage_candidate` (task 23 §5.2 — flagging, not filtering).
    emit == 'single': the single best by (coverage, identity, frac).
    Non-selected candidates are reclassified `secondary_target`.
    """
    candidates = [
        r for r in rows if r['classification'] == 'target_candidate']
    non_candidates = [
        r for r in rows if r['classification'] != 'target_candidate']

    if not candidates:
        return [], rows

    ranked = sorted(
        candidates,
        key=lambda r: (
            r['coverage'], r['ref_identity_pct'], r['aligned_frac']),
        reverse=True,
    )

    if thresholds['emit'] == 'all':
        keep = ranked[:thresholds['max_contigs']]
        floor = thresholds['low_coverage_fraction'] * keep[0]['coverage']
        primaries = [
            {**r, 'low_coverage_candidate': r['coverage'] < floor}
            for r in keep
        ]
        dropped = [
            {**r, 'classification': 'secondary_target'}
            for r in ranked[thresholds['max_contigs']:]
        ]
        return primaries, non_candidates + dropped

    primaries = ranked[:1]
    runners_up = [
        {**r, 'classification': 'secondary_target'} for r in ranked[1:]
    ]
    return primaries, non_candidates + runners_up


def write_secondaries(path: Path, rows: list, panels: list) -> None:
    """Write secondaries.tsv with one aligned-fraction column per panel."""
    frac_cols = [f'aligned_frac_{panel}' for panel in panels]
    header = BASE_SECONDARIES_COLS[:-1] + frac_cols + ['classification']
    with open(path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for row in rows:
            values = {
                **row,
                **{
                    f'aligned_frac_{panel}':
                        row['panel_aligned_frac'].get(panel, '')
                    for panel in panels
                },
            }
            fh.write('\t'.join(str(values[c]) for c in header) + '\n')


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assembly', type=Path, required=True)
    p.add_argument('--assembly-info', type=Path, required=True)
    p.add_argument('--gfa', type=Path, required=True)
    p.add_argument('--ref-dir', type=Path, required=True,
                   help='Reference bundle directory holding ${target}.mmi')
    p.add_argument('--sample-id', required=True)
    p.add_argument('--assembly-target', required=True)
    p.add_argument('--genetic-codes', required=True,
                   help='Comma-separated NCBI genetic-code table IDs')
    p.add_argument('--min-identity', type=float, required=True,
                   help='Merged identity floor, percent (spec §3.7.6)')
    p.add_argument('--min-aligned-frac', type=float, required=True,
                   help='Merged aligned-fraction floor (spec §3.7.6)')
    p.add_argument('--emit', choices=['single', 'all'], required=True,
                   help='Emit the top candidate, or all of them')
    p.add_argument('--max-contigs', type=int, required=True,
                   help='Cap on emitted contigs when --emit all')
    p.add_argument('--low-coverage-fraction', type=float, required=True,
                   help='Flag candidates below this fraction of the '
                        'top candidate coverage (task 23 §5.2)')
    p.add_argument('--out-target', type=Path, required=True)
    p.add_argument('--out-secondaries', type=Path, required=True)
    p.add_argument('--out-metadata', type=Path, required=True)
    args = p.parse_args()

    tables = [int(t) for t in args.genetic_codes.split(',')]
    thresholds = {
        'min_identity': args.min_identity,
        'min_aligned_frac': args.min_aligned_frac,
        'emit': args.emit,
        'max_contigs': args.max_contigs,
        'low_coverage_fraction': args.low_coverage_fraction,
    }

    info = parse_assembly_info(args.assembly_info)
    sequences = {
        r.id: str(r.seq)
        for r in SeqIO.parse(args.assembly, 'fasta')
    }

    aligners, sibling_panels_scored = load_panels(
        args.ref_dir, args.assembly_target)
    panels_scored = [args.assembly_target] + sibling_panels_scored

    rows = classify_contigs(
        info, sequences, aligners, args.assembly_target, tables,
        thresholds)
    primaries, secondaries = select_primary(
        rows, args.assembly_target, thresholds)

    with open(args.out_target, 'w') as fh:
        for row in primaries:
            cid = row['contig_id']
            seq = sequences.get(cid, '')
            fh.write(f'>{cid}\n{seq}\n')

    write_secondaries(args.out_secondaries, secondaries, panels_scored)

    audit_cols = [
        'contig_id', 'length_bp', 'coverage', 'ref_identity_pct',
        'aligned_frac', 'panel_aligned_frac', 'panel_identity_pct',
        'winning_panel', 'circular', 'circular_method', 'orf_aa_len',
        'selected_genetic_code', 'low_coverage_candidate',
        'classification',
    ]
    metadata = {
        'sample_id': args.sample_id,
        'assembly_target': args.assembly_target,
        'thresholds_applied': thresholds,
        'panels_scored': panels_scored,
        'sibling_panels_scored': sibling_panels_scored,
        'n_contigs_total': len(rows),
        'n_target_selected': len(primaries),
        'contigs_selected': [r['contig_id'] for r in primaries],
        'contigs': [
            {c: row[c] for c in audit_cols}
            for row in primaries + secondaries
        ],
    }
    if primaries:
        best = primaries[0]
        metadata['selected_genetic_code'] = best['selected_genetic_code']
        metadata['ref_identity_pct'] = best['ref_identity_pct']
        metadata['coverage'] = best['coverage']
        metadata['orf_aa_len'] = best['orf_aa_len']
        metadata['circular'] = best['circular']
        metadata['circular_method'] = best['circular_method']

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
