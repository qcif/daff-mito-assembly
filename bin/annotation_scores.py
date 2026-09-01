#!/usr/bin/env python3
"""SCORE_ANNOTATION step 3/3 — join blastp scores back on (task 40).

Reads the blastp TSV (``qseqid sseqid pident qcovhsp evalue
bitscore``) produced against the same per-gene protein panel
``MINIPROT_CDS`` aligned to, picks the best hit per feature (highest
bitscore; ties broken on identity), and attaches ``pident`` /
``qcovhsp`` / ``bitscore`` to every CDS feature in the merged GFF and
to ``annotation_summary.json``'s new ``cds_scores`` list — whichever
tool called the feature, and whether or not it has a hit at all.

Scores are data, never a gate (task 40 §5): no feature is ever
dropped, re-called, or excluded from the GFF/summary here, and this
module never reads or reasons about the summary's ``status`` field.
A feature with no blastp hit gets an explicit ``null`` triplet, not a
zero — a missing hit and a zero-identity hit are not the same thing,
and collapsing them would hide the difference from a reviewer.

Always exits 0 (CONSTITUTION rule 8).
"""

import argparse
import json
import sys
from pathlib import Path

from annotation_gff import parse_annotation_gff, parse_attrs

SCHEMA = "wf5/annotation-summary/v3"

BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "qcovhsp", "evalue", "bitscore",
]
SCORE_FIELDS = ("pident", "qcovhsp", "bitscore")

NULL_SCORE = {"pident": None, "qcovhsp": None, "bitscore": None,
              "sseqid": None}


def parse_blast_tsv(path: Path) -> dict:
    """Best hit per qseqid: highest bitscore, ties broken on pident.
    Malformed/truncated rows are skipped with a warning, not fatal."""
    best = {}
    if not path.exists():
        return best
    with open(path) as fh:
        for line_no, raw_line in enumerate(fh, 1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != len(BLAST_COLUMNS):
                print(
                    f"WARNING: skipping malformed blast line {line_no} "
                    f"in {path}: expected {len(BLAST_COLUMNS)} columns, "
                    f"got {len(fields)}",
                    file=sys.stderr,
                )
                continue
            qseqid, sseqid, pident, qcovhsp, _evalue, bitscore = fields
            try:
                pident = float(pident)
                qcovhsp = float(qcovhsp)
                bitscore = float(bitscore)
            except ValueError as exc:
                print(
                    f"WARNING: skipping malformed blast line {line_no} "
                    f"in {path}: {exc}",
                    file=sys.stderr,
                )
                continue
            current = best.get(qseqid)
            if current is None or (bitscore, pident) > (
                    current["bitscore"], current["pident"]):
                best[qseqid] = {
                    "pident": pident, "qcovhsp": qcovhsp,
                    "bitscore": bitscore, "sseqid": sseqid,
                }
    return best


def format_attr_value(value) -> str:
    return "null" if value is None else str(value)


def rescored_line(line: str, scores_by_id: dict) -> str:
    """Return `line` unchanged unless it's an mRNA row for a scored
    feature, in which case the score triplet is appended verbatim to
    its attributes."""
    fields = line.split("\t")
    if len(fields) != 9 or fields[2] != "mRNA":
        return line
    fid = parse_attrs(fields[8]).get("ID")
    if fid is None:
        return line
    score = scores_by_id.get(fid, NULL_SCORE)
    suffix = ";".join(
        f"{field}={format_attr_value(score[field])}"
        for field in SCORE_FIELDS
    )
    fields[8] = fields[8].rstrip(";") + ";" + suffix
    return "\t".join(fields)


def write_scored_gff(
    merged_gff: Path, out_gff: Path, scores_by_id: dict,
) -> None:
    lines = merged_gff.read_text().split("\n")
    scored = [rescored_line(line, scores_by_id) for line in lines]
    out_gff.write_text("\n".join(scored))


def cds_scores(features: list, scores_by_id: dict) -> list:
    rows = []
    for feature in sorted(
            features, key=lambda f: (f["seqid"], f["start"])):
        score = scores_by_id.get(feature["id"], NULL_SCORE)
        rows.append({
            "id": feature["id"],
            "gene": feature["gene"],
            "seqid": feature["seqid"],
            "start": feature["start"],
            "end": feature["end"],
            "source": feature["source"],
            "pident": score["pident"],
            "qcovhsp": score["qcovhsp"],
            "bitscore": score["bitscore"],
            "sseqid": score["sseqid"],
        })
    return rows


def run(
    merged_gff: Path, annotation_summary: Path, blast_tsv: Path,
    out_gff: Path, out_summary: Path,
) -> None:
    features = parse_annotation_gff(merged_gff)
    scores_by_id = parse_blast_tsv(blast_tsv)

    write_scored_gff(merged_gff, out_gff, scores_by_id)

    summary = json.loads(annotation_summary.read_text())
    summary["$schema"] = SCHEMA
    summary["cds_scores"] = cds_scores(features, scores_by_id)
    out_summary.write_text(json.dumps(summary, indent=2) + "\n")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--merged-gff", type=Path, required=True)
    p.add_argument("--annotation-summary", type=Path, required=True)
    p.add_argument("--blast-tsv", type=Path, required=True)
    p.add_argument("--out-gff", type=Path, required=True)
    p.add_argument("--out-summary", type=Path, required=True)
    args = p.parse_args()

    run(
        args.merged_gff, args.annotation_summary, args.blast_tsv,
        args.out_gff, args.out_summary,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
