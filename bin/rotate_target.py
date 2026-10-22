#!/usr/bin/env python3
"""ROTATE_TARGET — spec §2 stage 9a, §2.2 C13, §3.8, task 76.

Flye writes a circular mitogenome as a linear record cut at an arbitrary
origin that can differ between runs. This component puts each circular
`animal_mt` target contig into a fixed frame: it finds the first anchor
gene (COX1, then the fallbacks) from ROTATE_ANCHORS' miniprot output,
reverse-complements the contig when the anchor is on the minus strand,
and rotates it so the anchor CDS starts at position 1.

Rotation is a permutation (reverse-complement and cyclic shift). When
tandem parameters are given, a circular contig that holds the genome
two or more times in a row is then cut back to one copy, starting at the
anchor (task 75). FASTA header lines are written back exactly as read;
provenance goes to ``rotation.json`` (schema ``wf5/rotation/v2``).
Anything short of unreadable input is a passthrough with a recorded
status (CONSTITUTION principle 8).
"""

import argparse
import json
import sys
from pathlib import Path

SCHEMA = "wf5/rotation/v2"

STATUS_ROTATED = "rotated"
STATUS_NOT_CONFIGURED = "not_configured"
STATUS_NO_TARGET = "no_target"
STATUS_NOT_CIRCULAR = "not_circular"
STATUS_NO_ANCHOR = "no_anchor"

TANDEM_NOT_TESTED = "not_tested"
TANDEM_SINGLE = "single"
TANDEM_COLLAPSED = "collapsed"
TANDEM_UNRESOLVED = "unresolved"

CUT_LEADING_KMER = "leading_kmer"
CUT_PERIOD = "period"

# Upper bound on how many times a k-mer may recur and still be indexed
# for the period search.
MAX_COPIES = 4
TANDEM_PARAM_KEYS = ("kmer", "min_support", "tolerance_bp", "min_period_bp")

CIRCULAR_METHODS = ("flye_circ", "end_overlap")
CODON_LEN = 3
LINE_WIDTH = 80

COMPLEMENT = str.maketrans(
    "ACGTUNRYKMSWBDHVacgtunrykmswbdhv",
    "TGCAANYRMKSWVHDBtgcaanyrmkswvhdb",
)


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path: Path) -> list:
    """Return [(header_line, sequence), ...]; header keeps its '>' and
    any description exactly as written."""
    records = []
    header = None
    chunks = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def contig_id_of(header: str) -> str:
    return header[1:].split()[0] if len(header) > 1 else ""


def write_fasta(path: Path, records: list) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(header + "\n")
            for i in range(0, len(seq), LINE_WIDTH):
                fh.write(seq[i:i + LINE_WIDTH] + "\n")


def parse_attrs(field: str) -> dict:
    return dict(
        kv.split("=", 1) for kv in field.strip().split(";") if "=" in kv
    )


def parse_hits(path: Path, table: int, anchors: list) -> list:
    """Return anchor-gene hits from one miniprot GFF3's mRNA rows.

    Each hit is a dict with contig, gene, table, score, start, end
    (1-based, inclusive), strand and query_start. Unusable rows are
    skipped."""
    hits = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9 or fields[2] != "mRNA":
                continue
            target = parse_attrs(fields[8]).get("Target", "").split()
            if len(target) < 2:
                continue
            gene = target[0].rsplit("_", 1)[-1]
            if gene not in anchors or fields[6] not in "+-":
                continue
            try:
                hit = {
                    "contig": fields[0],
                    "gene": gene,
                    "table": table,
                    "score": float(fields[5]),
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "query_start": int(target[1]),
                }
            except ValueError:
                continue
            hits.append(hit)
    return hits


def best_hit(hits: list):
    """Highest score, then lowest start, then lowest table."""
    return min(
        hits,
        key=lambda h: (-h["score"], h["start"], h["table"]),
        default=None,
    )


def is_circular(contig_meta: dict) -> bool:
    return (
        bool(contig_meta.get("circular"))
        and contig_meta.get("circular_method") in CIRCULAR_METHODS
    )


def choose_anchor(contig_hits: list, anchors: list):
    """First gene in anchor order that has any hit; order beats score."""
    for gene in anchors:
        hit = best_hit([h for h in contig_hits if h["gene"] == gene])
        if hit:
            return hit
    return None


def rotate_contig(seq: str, hit: dict) -> tuple:
    """Return (rotated_seq, details) for the winning anchor hit.

    The start estimate is the first base of the anchor CDS on its strand,
    extended upstream by 3 x (query_start - 1) bp and wrapped around the
    circle. ``raw_start`` is the 1-based first aligned base in the
    orientation that was applied."""
    length = len(seq)
    reverse = hit["strand"] == "-"
    if reverse:
        seq = revcomp(seq)
        raw_start = length - hit["end"] + 1
    else:
        raw_start = hit["start"]
    extension = CODON_LEN * (hit["query_start"] - 1)
    offset = (raw_start - 1 - extension) % length
    rotated = seq[offset:] + seq[:offset]
    return rotated, {
        "reverse_complemented": reverse,
        "raw_start": raw_start,
        "query_start": hit["query_start"],
        "offset": offset,
    }


def circular_distance(a: int, b: int, length: int) -> int:
    delta = abs(a - b)
    return min(delta, length - delta)


def index_circular_kmers(seq: str, k: int) -> dict:
    """Map each k-mer, read around the circle, to its start positions;
    only k-mers that occur 2..MAX_COPIES times are kept."""
    wrapped = seq + seq[:k - 1]
    positions = {}
    for i in range(len(seq)):
        positions.setdefault(wrapped[i:i + k], []).append(i)
    return {
        kmer: pos for kmer, pos in positions.items()
        if 2 <= len(pos) <= MAX_COPIES
    }


def find_period(seq: str, k: int, tolerance_bp: int) -> tuple:
    """Return (period, support, copies) for a circular sequence.

    The period is the most common circular distance between adjacent
    occurrences of a repeated k-mer; support is the fraction of positions
    whose k-mer recurs at that distance (+/- tolerance_bp). Returns
    (0, 0.0, 0) when no k-mer repeats."""
    length = len(seq)
    repeated = index_circular_kmers(seq, k)
    histogram = {}
    for pos in repeated.values():
        for a, b in zip(pos, pos[1:] + pos[:1]):
            d = circular_distance(a, b, length)
            histogram[d] = histogram.get(d, 0) + 1
    if not histogram:
        return 0, 0.0, 0
    period = max(histogram, key=lambda d: (histogram[d], -d))
    supported = 0
    for pos in repeated.values():
        for a in pos:
            if any(
                abs(circular_distance(a, b, length) - period)
                <= tolerance_bp
                for b in pos if b != a
            ):
                supported += 1
    return period, supported / length, round(length / period)


def tandem_decision(
    length: int, period: int, support: float, copies: int, params: dict,
) -> str:
    """Apply the task 75 decision table to a measured period."""
    if (
        period == 0
        or support < params["min_support"]
        or period < params["min_period_bp"]
    ):
        return TANDEM_SINGLE
    if (
        copies >= 2
        and abs(length - copies * period)
        <= copies * params["tolerance_bp"]
    ):
        return TANDEM_COLLAPSED
    return TANDEM_UNRESOLVED


def find_cut(rotated: str, period: int, k: int, tolerance_bp: int) -> tuple:
    """Return (cut, method): where the first copy ends in the rotated
    frame. The cut is the next occurrence of the leading k-mer near the
    period, else the period itself."""
    lead = rotated[:k]
    cut = rotated.find(lead, max(period - tolerance_bp, 1))
    if cut != -1 and cut <= period + tolerance_bp:
        return cut, CUT_LEADING_KMER
    return period, CUT_PERIOD


def blank_tandem(status: str, params) -> dict:
    return {
        "status": status,
        "period_bp": None,
        "copies": None,
        "support": None,
        "cut_bp": None,
        "cut_method": None,
        "params": dict(params) if params else None,
    }


def detect_tandem(seq: str, circular: bool, params) -> dict:
    """Measure and classify one contig; the cut is added later."""
    if not params or not circular or not seq:
        return blank_tandem(TANDEM_NOT_TESTED, params)
    period, support, copies = find_period(
        seq, params["kmer"], params["tolerance_bp"])
    tandem = blank_tandem(
        tandem_decision(len(seq), period, support, copies, params), params)
    tandem.update({
        "period_bp": period or None,
        "copies": copies or None,
        "support": round(support, 4),
    })
    return tandem


def blank_record(contig_id: str, length: int, anchors: list) -> dict:
    return {
        "contig_id": contig_id,
        "status": None,
        "anchor_gene": None,
        "fallback": False,
        "anchor_table": None,
        "reverse_complemented": False,
        "raw_start": None,
        "query_start": None,
        "offset": 0,
        "length_bp": length,
        "input_length_bp": length,
        "tandem": None,
        "anchors": list(anchors),
    }


def rotate_record(
    record: dict, seq: str, circular: bool, hits: list, anchors: list,
) -> str:
    """Fill the rotation fields of ``record``; return the rotated seq."""
    if not anchors:
        record["status"] = STATUS_NOT_CONFIGURED
        return seq
    if not circular:
        record["status"] = STATUS_NOT_CIRCULAR
        return seq
    hit = choose_anchor(hits, anchors)
    if hit is None or not seq:
        record["status"] = STATUS_NO_ANCHOR
        return seq
    rotated, details = rotate_contig(seq, hit)
    record.update(details)
    record.update({
        "status": STATUS_ROTATED,
        "anchor_gene": hit["gene"],
        "fallback": hit["gene"] != anchors[0],
        "anchor_table": hit["table"],
    })
    return rotated


def process_contig(
    header: str, seq: str, contig_meta: dict, hits: list, anchors: list,
    tandem_params: dict = None,
) -> tuple:
    """Return (record, output_sequence) for one contig.

    Order: detect on the input, rotate, then cut to one copy so the kept
    copy starts at the anchor."""
    record = blank_record(contig_id_of(header), len(seq), anchors)
    circular = is_circular(contig_meta)
    tandem = detect_tandem(seq, circular, tandem_params)
    out = rotate_record(record, seq, circular, hits, anchors)
    if tandem["status"] == TANDEM_COLLAPSED:
        cut, method = find_cut(
            out, tandem["period_bp"], tandem_params["kmer"],
            tandem_params["tolerance_bp"])
        tandem.update({"cut_bp": cut, "cut_method": method})
        out = out[:cut]
        record["length_bp"] = cut
    record["tandem"] = tandem
    return record, out


def rotate_all(
    records: list, bin_metadata: dict, hits: list, anchors: list,
    tandem_params: dict = None,
) -> tuple:
    """Return (output_records, rotation_json)."""
    meta_by_id = {
        c.get("contig_id"): c for c in bin_metadata.get("contigs") or []
    }
    hits_by_contig = {}
    for hit in hits:
        hits_by_contig.setdefault(hit["contig"], []).append(hit)

    out_records = []
    contig_records = []
    for header, seq in records:
        cid = contig_id_of(header)
        record, new_seq = process_contig(
            header, seq, meta_by_id.get(cid, {}),
            hits_by_contig.get(cid, []), anchors, tandem_params)
        out_records.append((header, new_seq))
        contig_records.append(record)

    if not anchors:
        status = STATUS_NOT_CONFIGURED
    elif not records:
        status = STATUS_NO_TARGET
    else:
        status = max(
            zip(contig_records, records),
            key=lambda pair: len(pair[1][1]),
        )[0]["status"]
    return out_records, {
        "$schema": SCHEMA,
        "status": status,
        "anchors": list(anchors),
        "tandem_params": dict(tandem_params) if tandem_params else None,
        "contigs": contig_records,
    }


def parse_gff_arg(raw: str) -> tuple:
    table_str, path_str = raw.split(":", 1)
    return int(table_str), Path(path_str)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target-fasta", type=Path, required=True)
    p.add_argument("--bin-metadata", type=Path, required=True)
    p.add_argument(
        "--anchors", default="",
        help="comma-separated anchor genes in priority order; empty "
             "disables rotation")
    p.add_argument(
        "--gff", action="append", default=[], metavar="TABLE:PATH",
        help="ROTATE_ANCHORS miniprot GFF for one genetic-code table")
    p.add_argument(
        "--tandem-kmer", type=int,
        help="k-mer size for tandem-duplicate detection; omit all four "
             "--tandem-* args to disable detection")
    p.add_argument("--tandem-min-support", type=float)
    p.add_argument("--tandem-tolerance-bp", type=int)
    p.add_argument("--tandem-min-period-bp", type=int)
    p.add_argument("--out-fasta", type=Path, required=True)
    p.add_argument("--out-json", type=Path, required=True)
    args = p.parse_args()

    anchors = [g for g in args.anchors.split(",") if g]
    records = read_fasta(args.target_fasta)
    bin_metadata = json.loads(args.bin_metadata.read_text())

    hits = []
    for raw in args.gff:
        table, path = parse_gff_arg(raw)
        hits.extend(parse_hits(path, table, anchors))

    tandem_values = [
        args.tandem_kmer, args.tandem_min_support,
        args.tandem_tolerance_bp, args.tandem_min_period_bp,
    ]
    tandem_params = None
    if all(v is not None for v in tandem_values):
        tandem_params = dict(zip(TANDEM_PARAM_KEYS, tandem_values))

    out_records, rotation = rotate_all(
        records, bin_metadata, hits, anchors, tandem_params)
    write_fasta(args.out_fasta, out_records)
    args.out_json.write_text(json.dumps(rotation, indent=2) + "\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
