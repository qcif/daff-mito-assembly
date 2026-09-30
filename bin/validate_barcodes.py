#!/usr/bin/env python3
"""EXTRACT_BARCODES (C5) — spec §2 stage 13, §2.2 C5, §3.3.

Subsets ``MINIPROT_CDS``'s ``cds.gff`` to the loci named in the locus
panel (``assets/loci.json``), then validates each candidate: length,
ORF integrity under the target-appropriate NCBI genetic-code table
(clade trial on ``animal_mt``), internal stop codons, and a
protein-identity floor.

Never re-aligns: every emitted barcode's coordinates are either copied
verbatim from its ``cds.gff`` source row(s) (``pass``), or a
verbatim, contiguous sub-span of that row shipped as a labelled
``partial`` (task 51 — a strongly-supported locus that fails only on a
short, isolated internal stop is trimmed to its longest stop-free
stretch rather than dropped outright, per CONSTITUTION principle 7).
That trim only ever narrows the emitted span; it never changes a base,
never touches a passing locus, and is recorded in full in
``validation.tsv`` and ``coords.gff`` (principle 18). The ORF check
otherwise never changes what gets extracted or emitted, only how it is
validated.

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
nucleotide sequence is untouched either way; only the pass/partial/fail
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

STATUS_PASS = 'pass'
STATUS_PARTIAL = 'partial'
STATUS_FAIL = 'fail'
STATUS_NOT_FOUND = 'not_found'

# Taxodactyl's own accepted-partial floor (confirmed with the user,
# 2026-09-29) — the *consumer's* hard minimum, not a target. This
# workflow's own --partial-min-nt default (nextflow.config) is set well
# above it on purpose, so this constant only ever rejects a config
# value that would ship a barcode the consumer will refuse (task 51
# §3.4, rule 18).
TAXODACTYL_MIN_PARTIAL_NT = 20

TSV_COLUMNS = [
    'gene', 'status', 'reason', 'seqid', 'start', 'end', 'strand',
    'identity', 'genetic_code', 'length_nt',
    'source_start', 'source_end', 'source_length_nt', 'n_internal_stops',
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
                        'mrna_id': attrs['ID'],
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


def codon_blocks_with_offsets(seq: str, cigar: str) -> list:
    """Like ``codon_blocks``, but each block also carries the offset
    (into ``seq``) of its first base — needed to map a trimmed segment
    (§3.2) back to genome coordinates. The two functions share one
    walk so there is exactly one definition of a codon boundary
    (CONSTITUTION rule 19); ``codon_blocks`` is a thin view over this.
    """
    blocks = []
    current = ''
    current_offset = None
    cursor = 0
    for n_str, op in CIGAR_OP_RE.findall(cigar):
        n = int(n_str)
        if op in FRAME_SAFE_OPS:
            if current_offset is None:
                current_offset = cursor
            current += seq[cursor:cursor + n * 3]
            cursor += n * 3
        elif op == 'I':
            continue
        else:  # BLOCK_BREAK_OPS: frameshift or intron — drop, break
            if current:
                blocks.append((current_offset, current))
                current = ''
                current_offset = None
            cursor += n
    remainder = seq[cursor:]
    if remainder and len(remainder) % 3 == 0:
        if current_offset is None:
            current_offset = cursor
        current += remainder
    if current:
        blocks.append((current_offset, current))
    return blocks


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
    return [block for _offset, block in codon_blocks_with_offsets(seq, cigar)]


def build_blocks(seq: str, cigar) -> list:
    """Codon blocks with offsets (§3.2), from a CIGAR when there is
    one, otherwise the whole sequence truncated to a full number of
    codons starting at offset 0."""
    if cigar:
        return codon_blocks_with_offsets(seq, cigar)
    trunc = seq[:len(seq) - len(seq) % 3]
    return [(0, trunc)] if trunc else []


def block_has_internal_stop(block: str, table: int, is_last: bool) -> bool:
    protein = str(Seq(block).translate(table=table, to_stop=False))
    if is_last and protein.endswith('*'):
        protein = protein[:-1]
    return '*' in protein


def find_orf_segments(blocks: list, table: int):
    """Split the translation of ``blocks`` under ``table`` at every
    internal stop codon. Frameshift/intron breaks between blocks never
    end a segment (§3.2) — only a stop codon does, and its own three
    bases are excluded from both neighbouring segments. The normal
    terminal-stop allowance still applies: a stop as the final codon of
    the final block is not internal and stays inside the segment that
    reaches it.

    Returns ``(segments, n_internal_stops)``, where ``segments`` is a
    list of ``(start, end)`` exclusive-end spans in ``seq``-offset
    space, in 5'->3' order — so a tie-break on segment length that
    prefers the first one encountered is automatically the 5'-most
    (§3.2).
    """
    if not blocks:
        return [], 0
    segments = []
    n_stops = 0
    n_blocks = len(blocks)
    seg_start = blocks[0][0]
    for bi, (offset, block) in enumerate(blocks):
        is_last_block = bi == n_blocks - 1
        protein = str(Seq(block).translate(table=table, to_stop=False))
        n_codons = len(protein)
        for ci, aa in enumerate(protein):
            if aa != '*':
                continue
            is_terminal = is_last_block and ci == n_codons - 1
            if is_terminal:
                continue
            n_stops += 1
            codon_start = offset + ci * 3
            segments.append((seg_start, codon_start))
            seg_start = codon_start + 3
    # The terminal codon of the last block is never itself treated as
    # a stop-triggering reset point (see the ``continue`` above), so
    # at least it always remains open past the last reset — this final
    # segment is never empty.
    block_end = blocks[-1][0] + len(blocks[-1][1])
    segments.append((seg_start, block_end))
    return segments, n_stops


def spliced_exon_offsets(feature: dict) -> list:
    """The same exon ordering ``extract_cds_seq`` splices in (ascending
    start, reversed on minus strand), annotated with each exon's
    ``seq``-offset span alongside its own genome coordinates."""
    ordered = sorted(
        feature['cds'], key=lambda iv: iv[0],
        reverse=(feature['strand'] == '-'),
    )
    offsets = []
    cursor = 0
    for start, end in ordered:
        length = end - start + 1
        offsets.append((cursor, cursor + length, start, end))
        cursor += length
    return offsets


def map_segment_to_genome(feature: dict, seg_start: int, seg_end: int):
    """Map a ``[seg_start, seg_end)`` span of the spliced CDS sequence
    back to 1-based inclusive genome coordinates. Returns ``None`` if
    the span would cross a CDS exon junction — a multi-exon partial is
    not built (§3.2); the caller fails the locus instead of guessing a
    spliced sub-range."""
    for off_start, off_end, exon_start, exon_end in \
            spliced_exon_offsets(feature):
        if seg_start >= off_start and seg_end <= off_end:
            local_start = seg_start - off_start
            local_end = seg_end - off_start
            if feature['strand'] == '-':
                return exon_end - local_end + 1, exon_end - local_start
            return exon_start + local_start, exon_start + local_end - 1
    return None


def validate_orf(
    seq, cigar, tables: list, min_identity: float, identity: float,
    feature: dict = None, partial_min_nt: int = 0,
):
    """Length, identity-floor and ORF/internal-stop checks, in that
    order (§3.1 — length and identity failures are never rescued).
    Returns a dict:

        status              'pass' | 'partial' | 'fail'
        reason               None, or one of the REASON_* constants
        table                chosen genetic-code table, or None
        n_internal_stops     int, or None where the TSV leaves it blank
        segment              (start, end) emitted span in seq-offset
                              space (whole sequence for 'pass')
        genome_span          (start, end) 1-based inclusive emitted
                              genome coordinates for 'partial'; None
                              otherwise (pass/fail reuse the feature's
                              own start/end untouched)

    With a usable CIGAR, translation is segmented at frameshift/intron
    breakpoints (find_orf_segments) so a single-base indel — the
    dominant real-world failure mode on ONT-derived assemblies (task 30
    outcomes) — doesn't corrupt the reading frame for the rest of the
    gene. Without one (e.g. a hand-built cds.gff with no ##PAF
    context), falls back to translating the whole sequence truncated
    to a whole number of codons.

    When every configured table has an internal stop and
    ``partial_min_nt`` is nonzero, the table+segment giving the longest
    stop-free stretch is used to ship a labelled ``partial`` instead of
    failing outright, provided that stretch is at least
    ``partial_min_nt`` and maps cleanly to a single CDS exon (task 51).
    """
    empty_result = {
        'status': STATUS_FAIL, 'reason': None, 'table': None,
        'n_internal_stops': None, 'segment': None, 'genome_span': None,
    }
    if len(seq) == 0:
        return {**empty_result, 'reason': REASON_INVALID_LENGTH}
    if identity * 100 < min_identity:
        return {**empty_result, 'reason': REASON_IDENTITY_BELOW_FLOOR}

    blocks = build_blocks(str(seq), cigar)
    if not blocks:
        return {**empty_result, 'reason': REASON_INVALID_LENGTH}

    for table in tables:
        _segments, n_stops = find_orf_segments(blocks, table)
        if n_stops == 0:
            return {
                'status': STATUS_PASS, 'reason': None, 'table': table,
                'n_internal_stops': 0, 'segment': (0, len(seq)),
                'genome_span': None,
            }

    # Every table has at least one internal stop. n_internal_stops on
    # a plain fail is always counted under the first configured table
    # (§3.4), independent of which table the rescue attempt below ends
    # up preferring.
    _, first_table_stops = find_orf_segments(blocks, tables[0])

    if partial_min_nt > 0 and feature is not None:
        # tables (a required CLI arg) and blocks (checked above) are
        # both non-empty, and find_orf_segments always yields at least
        # one segment for non-empty blocks, so a best candidate always
        # exists after this loop.
        best = None  # (table, (start, end), n_stops)
        for table in tables:
            segments, n_stops = find_orf_segments(blocks, table)
            candidate = max(segments, key=lambda s: s[1] - s[0])
            candidate_len = candidate[1] - candidate[0]
            if best is None or candidate_len > best[1][1] - best[1][0]:
                best = (table, candidate, n_stops)
        table, segment, n_stops = best
        if segment[1] - segment[0] >= partial_min_nt:
            genome_span = map_segment_to_genome(feature, *segment)
            # genome_span is None only if the segment would cross a
            # CDS exon junction — don't guess a spliced sub-range,
            # fail below instead (§3.2).
            if genome_span is not None:
                return {
                    'status': STATUS_PARTIAL,
                    'reason': REASON_INTERNAL_STOP,
                    'table': table, 'n_internal_stops': n_stops,
                    'segment': segment, 'genome_span': genome_span,
                }

    return {
        **empty_result,
        'reason': REASON_INTERNAL_STOP,
        'n_internal_stops': first_table_stops,
    }


def tsv_row(
    gene, status, reason='', feature=None, table=None, length_nt='',
    start=None, end=None, source_start='', source_end='',
    source_length_nt='', n_internal_stops='',
):
    row = {col: '' for col in TSV_COLUMNS}
    row['gene'] = gene
    row['status'] = status
    row['reason'] = reason
    row['genetic_code'] = table if table is not None else ''
    row['length_nt'] = length_nt
    row['source_start'] = source_start
    row['source_end'] = source_end
    row['source_length_nt'] = source_length_nt
    row['n_internal_stops'] = (
        n_internal_stops if n_internal_stops is not None else '')
    if feature is not None:
        row['seqid'] = feature['seqid']
        row['strand'] = feature['strand']
        row['identity'] = feature['identity']
        row['start'] = feature['start'] if start is None else start
        row['end'] = feature['end'] if end is None else end
    return row


def build_partial_coords_line(feature: dict, start: int, end: int) -> str:
    """One derived GFF line for a partial's emitted span, alongside the
    source cds.gff row(s) reproduced verbatim (§3.4) — links back to
    the source mRNA feature and marks itself as a trim, never silent
    (CONSTITUTION rule 18)."""
    attrs = (
        f"Parent={feature['mrna_id']};partial=true;"
        f"trim_reason={REASON_INTERNAL_STOP}"
    )
    return (
        f"{feature['seqid']}\tvalidate_barcodes\tCDS\t{start}\t{end}\t"
        f".\t{feature['strand']}\t.\t{attrs}"
    )


def run(cds_gff: Path, target_fasta: Path, assembly_target: str,
        locus_panel: Path, tables: list, min_identity: float,
        out_fasta: Path, out_coords: Path, out_tsv: Path,
        partial_min_nt: int = 0) -> None:
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
            result = validate_orf(
                seq, feature['cigar'], tables, min_identity,
                feature['identity'], feature=feature,
                partial_min_nt=partial_min_nt)
            status = result['status']

            if status == STATUS_FAIL:
                tsv_rows.append(tsv_row(
                    locus, 'fail', reason=result['reason'],
                    feature=feature, length_nt=len(seq),
                    n_internal_stops=result['n_internal_stops']))
                continue

            if status == STATUS_PASS:
                barcode_id = (
                    f"{locus}_{feature['seqid']}_"
                    f"{feature['start']}_{feature['end']}"
                )
                fasta_lines.append(f'>{barcode_id}')
                fasta_lines.append(str(seq))
                coords_lines.extend(feature['cds_lines'])
                tsv_rows.append(tsv_row(
                    locus, 'pass', feature=feature, table=result['table'],
                    length_nt=len(seq), source_start=feature['start'],
                    source_end=feature['end'], source_length_nt=len(seq),
                    n_internal_stops=0))
                continue

            # status == STATUS_PARTIAL
            seg_start, seg_end = result['segment']
            g_start, g_end = result['genome_span']
            emitted_seq = seq[seg_start:seg_end]
            barcode_id = (
                f"{locus}_{feature['seqid']}_{g_start}_{g_end}"
            )
            fasta_lines.append(
                f'>{barcode_id} partial=internal_stop_trimmed '
                f"source={feature['start']}-{feature['end']}")
            fasta_lines.append(str(emitted_seq))
            coords_lines.extend(feature['cds_lines'])
            coords_lines.append(
                build_partial_coords_line(feature, g_start, g_end))
            tsv_rows.append(tsv_row(
                locus, 'partial', reason=result['reason'],
                feature=feature, table=result['table'],
                length_nt=len(emitted_seq), start=g_start, end=g_end,
                source_start=feature['start'], source_end=feature['end'],
                source_length_nt=len(seq),
                n_internal_stops=result['n_internal_stops']))

    out_fasta.write_text(
        '\n'.join(fasta_lines) + ('\n' if fasta_lines else ''))
    out_coords.write_text(
        '\n'.join(coords_lines) + ('\n' if coords_lines else ''))

    with open(out_tsv, 'w', newline='') as fh:
        writer = csv.DictWriter(
            fh, fieldnames=TSV_COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(tsv_rows)


def partial_min_nt_type(value: str) -> int:
    """0 disables partial recovery entirely; any other value must clear
    Taxodactyl's own accepted-partial floor (task 51 §3.4) — shipping a
    partial the consumer will reject is worse than not shipping one
    (CONSTITUTION rule 18)."""
    n = int(value)
    if n != 0 and n < TAXODACTYL_MIN_PARTIAL_NT:
        raise argparse.ArgumentTypeError(
            '--partial-min-nt must be 0 (disabled) or at least '
            f'{TAXODACTYL_MIN_PARTIAL_NT} nt, the accepted-partial '
            f'floor Taxodactyl itself enforces; got {n}'
        )
    return n


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
    p.add_argument('--partial-min-nt', type=partial_min_nt_type,
                   default=0,
                   help='Minimum stop-free segment length (nt) to ship '
                        'as a labelled partial when an otherwise-good '
                        'locus fails only on an internal stop; 0 '
                        'disables the rescue (task 51 §3.4)')
    p.add_argument('--out-fasta', type=Path, required=True)
    p.add_argument('--out-coords', type=Path, required=True)
    p.add_argument('--out-tsv', type=Path, required=True)
    args = p.parse_args()

    tables = [int(t) for t in args.genetic_codes.split(',')]
    run(
        args.cds_gff, args.target_fasta, args.assembly_target,
        args.locus_panel, tables, args.min_identity,
        args.out_fasta, args.out_coords, args.out_tsv,
        partial_min_nt=args.partial_min_nt,
    )
    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
