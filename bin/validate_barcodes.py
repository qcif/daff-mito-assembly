#!/usr/bin/env python3
"""EXTRACT_BARCODES (C5) — spec §2 stage 13, §2.2 C5, §3.3.

Subsets ``MINIPROT_CDS``'s ``cds.gff`` to the loci named in the locus
panel (``assets/loci.json``), then validates each candidate: length,
ORF integrity under the target-appropriate NCBI genetic-code table
(clade trial on ``animal_mt``), internal stop codons, and a
protein-identity floor.

Never re-aligns or re-trims: every emitted barcode's coordinates are
copied verbatim from its ``cds.gff`` source row(s) — that invariant is
what the unified miniprot pass (task 30) buys, and it is what
``coords.gff`` exists to make auditable. The ORF check below never
changes what gets extracted or emitted, only how it is validated.

The comprehensive protein panel carries several representative protein
sequences per gene (task 29), so more than one ``cds.gff`` feature
commonly maps to the same genomic locus. Features are clustered by
gene symbol and genomic overlap; the highest-identity feature in each
cluster is the locus's candidate. Non-overlapping clusters for the
same gene are genuine distinct loci (e.g. a plastid gene inside the
inverted repeat) and are never collapsed.

**CIGAR-aware translation.** A real ONT-derived assembly routinely
carries a single-base indel error somewhere inside an otherwise
correct gene (miniprot flags this as ``Frameshift=N``) — this is the
dominant case in practice, not an edge case (task 30 outcomes). Plain
``--gff`` output already carries the extended CIGAR miniprot used for
each hit, in a ``##PAF`` comment line immediately preceding the hit's
``mRNA`` row (``cg:Z:`` tag; op semantics: miniprot(1) manpage —
``nM``/``nD`` advance n*3 genome nt in-frame, ``nI`` consumes protein
residues only, ``nF``/``nG``/``nN``/``nU``/``nV`` mark a frameshift or
intron of n genome nt that cannot be losslessly assigned to a codon).
Translation is segmented at those breakpoints so one bad base doesn't
corrupt the reading frame for the rest of the gene — the emitted
nucleotide sequence is untouched either way; only the pass/fail
verdict and chosen genetic-code table use the CIGAR. A candidate with
no usable CIGAR (e.g. a hand-built GFF without ``##PAF`` context)
falls back to whole-sequence translation.

Always exits 0 (CONSTITUTION rule 8) — a sample with no recoverable
barcode is a ``not_found`` result per locus, not a pipeline error.
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

# miniprot(1) manpage cg:Z op vocabulary. M/D advance the genome in
# whole codons (frame-safe); F/G/N/U/V advance it by an amount that
# cannot be assigned to a codon without fabricating bases, so they are
# treated as block breaks; I consumes no genome bases.
CIGAR_OP_RE = re.compile(r'(\d+)([MIDFGNUV])')
FRAME_SAFE_OPS = ('M', 'D')
BLOCK_BREAK_OPS = ('F', 'G', 'N', 'U', 'V')

REASON_NOT_FOUND = 'not_found'
REASON_INVALID_LENGTH = 'invalid_length'
REASON_IDENTITY_BELOW_FLOOR = 'identity_below_floor'
REASON_INTERNAL_STOP = 'internal_stop_codon'

TSV_COLUMNS = [
    'gene', 'status', 'reason', 'seqid', 'start', 'end', 'strand',
    'identity', 'genetic_code', 'length_nt',
]


def parse_attrs(field: str) -> dict:
    attrs = {}
    for kv in field.strip().split(';'):
        if not kv or '=' not in kv:
            continue
        key, value = kv.split('=', 1)
        attrs[key] = value
    return attrs


def cigar_tag(paf_line: str):
    """Pull the cg:Z: extended-CIGAR tag out of a ##PAF comment line,
    or None if the line carries no such tag."""
    for field in paf_line.split('\t'):
        if field.startswith('cg:Z:'):
            return field[len('cg:Z:'):]
    return None


def parse_gff(path: Path) -> list:
    """Parse mRNA + CDS rows into one record per miniprot hit.

    Each hit's ##PAF comment line (emitted by plain --gff, no --aln
    needed) immediately precedes its mRNA row and carries the cg:Z
    CIGAR for that hit — captured here so validate_orf can translate
    around frameshift breakpoints instead of naively end-to-end.

    Malformed lines are skipped with a warning; the rest of the file
    is still processed (spec principle 8 — never abort on bad input).
    """
    mrna = {}
    pending_cigar = None
    with open(path) as fh:
        for line_no, raw_line in enumerate(fh, 1):
            line = raw_line.rstrip('\n')
            if line.startswith('##PAF'):
                pending_cigar = cigar_tag(line)
                continue
            if not line or line.startswith('#'):
                continue
            try:
                fields = line.split('\t')
                seqid, _src, ftype, start, end, _score, strand, \
                    _phase, attr_field = fields
                start, end = int(start), int(end)
                attrs = parse_attrs(attr_field)
                if ftype == 'mRNA':
                    target = attrs['Target'].split()
                    protein_id = target[0]
                    gene = protein_id.rsplit('_', 1)[-1]
                    mrna[attrs['ID']] = {
                        'seqid': seqid,
                        'strand': strand,
                        'gene': gene,
                        'protein_id': protein_id,
                        'identity': float(attrs.get('Identity', 0.0)),
                        'cigar': pending_cigar,
                        'cds': [],
                        'cds_lines': [],
                    }
                    pending_cigar = None
                elif ftype == 'CDS':
                    parent = attrs['Parent']
                    if parent in mrna:
                        mrna[parent]['cds'].append((start, end))
                        mrna[parent]['cds_lines'].append(line)
            except (ValueError, KeyError, IndexError) as exc:
                print(
                    f'WARNING: skipping malformed GFF line {line_no}: '
                    f'{exc}',
                    file=sys.stderr,
                )
                continue
    return [rec for rec in mrna.values() if rec['cds']]


def cluster_by_locus(records: list) -> dict:
    """Group records by gene symbol, then merge genomically-overlapping
    hits (redundant panel representatives of the same locus), keeping
    the highest-identity hit per cluster. Non-overlapping hits are
    genuine distinct loci and are kept as separate clusters."""
    by_gene = {}
    for rec in records:
        rec['start'] = min(s for s, _e in rec['cds'])
        rec['end'] = max(e for _s, e in rec['cds'])
        by_gene.setdefault(rec['gene'].upper(), []).append(rec)

    winners_by_gene = {}
    for gene, recs in by_gene.items():
        recs = sorted(recs, key=lambda r: (r['seqid'], r['start']))
        clusters = []
        for rec in recs:
            if (clusters
                    and clusters[-1][-1]['seqid'] == rec['seqid']
                    and rec['start']
                    <= max(r['end'] for r in clusters[-1])):
                clusters[-1].append(rec)
            else:
                clusters.append([rec])
        winners_by_gene[gene] = [
            max(cluster, key=lambda r: r['identity'])
            for cluster in clusters
        ]
    return winners_by_gene


def extract_cds_seq(sequences: dict, feature: dict):
    """Splice the feature's CDS exons out of the target contig,
    respecting strand. Coordinates are 1-based inclusive (GFF)."""
    contig = sequences[feature['seqid']]
    ordered = sorted(
        feature['cds'], key=lambda iv: iv[0],
        reverse=(feature['strand'] == '-'),
    )
    parts = []
    for start, end in ordered:
        frag = contig[start - 1:end]
        if feature['strand'] == '-':
            frag = frag.reverse_complement()
        parts.append(frag)
    seq = parts[0]
    for part in parts[1:]:
        seq += part
    return seq


def codon_blocks(seq: str, cigar: str) -> list:
    """Split ``seq`` into frame-safe runs of whole codons, using the
    hit's cg:Z CIGAR to skip over frameshift/intron breakpoints that
    can't be assigned to a codon without inventing bases.

    ``seq`` is assumed to already be oriented 5'->3' (i.e. the CIGAR,
    which walks the alignment in that direction, can be applied
    front-to-back). Any trailing bases past the CIGAR's own span
    (typically a stop codon, appended to the alignment but not itself
    part of a protein residue) are folded into the last block when
    they form a whole codon.
    """
    blocks = []
    current = ''
    cursor = 0
    for n_str, op in CIGAR_OP_RE.findall(cigar):
        n = int(n_str)
        if op in FRAME_SAFE_OPS:
            current += seq[cursor:cursor + n * 3]
            cursor += n * 3
        elif op == 'I':
            continue
        else:  # BLOCK_BREAK_OPS: frameshift or intron — drop, break
            if current:
                blocks.append(current)
                current = ''
            cursor += n
    remainder = seq[cursor:]
    if remainder and len(remainder) % 3 == 0:
        current += remainder
    if current:
        blocks.append(current)
    return blocks


def block_has_internal_stop(block: str, table: int, is_last: bool) -> bool:
    protein = str(Seq(block).translate(table=table, to_stop=False))
    if is_last and protein.endswith('*'):
        protein = protein[:-1]
    return '*' in protein


def validate_orf(
    seq, cigar, tables: list, min_identity: float, identity: float,
):
    """Length, identity-floor and ORF/internal-stop checks, in that
    order. Returns (passed, reason, chosen_table).

    With a usable CIGAR, translation is segmented at frameshift/intron
    breakpoints (codon_blocks) so a single-base indel — the dominant
    real-world failure mode on ONT-derived assemblies (task 30
    outcomes) — doesn't corrupt the reading frame for the rest of the
    gene. Without one (e.g. a hand-built cds.gff with no ##PAF
    context), falls back to translating the whole sequence truncated
    to a whole number of codons.
    """
    if len(seq) == 0:
        return False, REASON_INVALID_LENGTH, None
    if identity * 100 < min_identity:
        return False, REASON_IDENTITY_BELOW_FLOOR, None
    if cigar:
        blocks = codon_blocks(str(seq), cigar)
    else:
        trunc = str(seq[:len(seq) - len(seq) % 3])
        blocks = [trunc] if trunc else []
    if not blocks:
        return False, REASON_INVALID_LENGTH, None
    for table in tables:
        if not any(
            block_has_internal_stop(b, table, i == len(blocks) - 1)
            for i, b in enumerate(blocks)
        ):
            return True, None, table
    return False, REASON_INTERNAL_STOP, None


def tsv_row(gene, status, reason='', feature=None, table=None,
            length_nt=''):
    row = {col: '' for col in TSV_COLUMNS}
    row['gene'] = gene
    row['status'] = status
    row['reason'] = reason
    row['genetic_code'] = table if table is not None else ''
    row['length_nt'] = length_nt
    if feature is not None:
        row['seqid'] = feature['seqid']
        row['start'] = feature['start']
        row['end'] = feature['end']
        row['strand'] = feature['strand']
        row['identity'] = feature['identity']
    return row


def run(cds_gff: Path, target_fasta: Path, assembly_target: str,
        locus_panel: Path, tables: list, min_identity: float,
        out_fasta: Path, out_coords: Path, out_tsv: Path) -> None:
    loci = json.loads(locus_panel.read_text())[assembly_target]

    sequences = {
        rec.id: rec.seq for rec in SeqIO.parse(target_fasta, 'fasta')
    }
    records = parse_gff(cds_gff)
    winners_by_gene = cluster_by_locus(records)

    tsv_rows = []
    fasta_lines = []
    coords_lines = []

    for locus in loci:
        winners = winners_by_gene.get(locus.upper(), [])
        if not winners:
            tsv_rows.append(
                tsv_row(locus, 'not_found', reason=REASON_NOT_FOUND))
            continue
        for feature in winners:
            seq = extract_cds_seq(sequences, feature)
            passed, reason, table = validate_orf(
                seq, feature['cigar'], tables, min_identity,
                feature['identity'])
            if not passed:
                tsv_rows.append(tsv_row(
                    locus, 'fail', reason=reason, feature=feature,
                    length_nt=len(seq)))
                continue
            barcode_id = (
                f"{locus}_{feature['seqid']}_"
                f"{feature['start']}_{feature['end']}"
            )
            fasta_lines.append(f'>{barcode_id}')
            fasta_lines.append(str(seq))
            coords_lines.extend(feature['cds_lines'])
            tsv_rows.append(tsv_row(
                locus, 'pass', feature=feature, table=table,
                length_nt=len(seq)))

    out_fasta.write_text(
        '\n'.join(fasta_lines) + ('\n' if fasta_lines else ''))
    out_coords.write_text(
        '\n'.join(coords_lines) + ('\n' if coords_lines else ''))

    with open(out_tsv, 'w', newline='') as fh:
        writer = csv.DictWriter(
            fh, fieldnames=TSV_COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(tsv_rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cds-gff', type=Path, required=True)
    p.add_argument('--target-fasta', type=Path, required=True)
    p.add_argument('--assembly-target', required=True)
    p.add_argument('--locus-panel', type=Path, required=True)
    p.add_argument('--genetic-codes', required=True,
                   help='Comma-separated NCBI genetic-code table IDs; '
                        'tried in order (clade trial on animal_mt)')
    p.add_argument('--min-identity', type=float, required=True,
                   help='Protein-identity floor, percent')
    p.add_argument('--out-fasta', type=Path, required=True)
    p.add_argument('--out-coords', type=Path, required=True)
    p.add_argument('--out-tsv', type=Path, required=True)
    args = p.parse_args()

    tables = [int(t) for t in args.genetic_codes.split(',')]
    run(
        args.cds_gff, args.target_fasta, args.assembly_target,
        args.locus_panel, tables, args.min_identity,
        args.out_fasta, args.out_coords, args.out_tsv,
    )
    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
