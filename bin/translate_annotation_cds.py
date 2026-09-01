#!/usr/bin/env python3
"""SCORE_ANNOTATION step 1/3 — translate CDS features (task 40 §2).

Reads the merged annotation GFF (``ANNOTATE``'s C8 output) and
translates every CDS feature — miniprot winner or MITOS2 rescue alike
— from its own genomic coordinates, under the sample's selected
genetic-code table (``genetic_code_cds`` in ``annotation_summary
.json``, task 38 §4). Both provenances go through exactly one
translation path, uniformly, from coordinates — miniprot's own
``--trans`` output is not used, since rescued features have no
equivalent and two translation paths would defeat the point of a
comparable score (task 40 §2).

Translation is whole-sequence, truncated to a full number of codons —
not segmented at frameshift/intron breakpoints the way
``validate_barcodes.codon_blocks`` does for barcode ORF validation.
That segmentation depends on the hit's own CIGAR, which a rescued
feature does not have; a heavily-frameshifted translation still
produces a garbled protein, but that is exactly what the following
blastp pass is for — it shows up as a poor identity/coverage score
rather than a silently "corrected" one (task 40 §1, ATP8 worked
example).

Runs in the wf5-scripts image (python3.12 + biopython) — the pinned
BLAST biocontainer that runs blastp next carries no Python at all
(task 40 §4).
"""

import argparse
import json
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from annotation_gff import parse_annotation_gff
from validate_barcodes import extract_cds_seq


def translate_feature(sequences: dict, feature: dict, table: int) -> str:
    """Splice + translate one feature's CDS, truncated to a whole
    number of codons. Returns '' if nothing is translatable (empty
    splice, or fewer than one full codon)."""
    seq = extract_cds_seq(sequences, feature)
    trunc = seq[:len(seq) - len(seq) % 3]
    if len(trunc) == 0:
        return ""
    protein = str(Seq(trunc).translate(table=table, to_stop=False))
    return protein.rstrip("*")


def run(
    merged_gff: Path, target_fasta: Path, annotation_summary: Path,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(annotation_summary.read_text())
    table = summary.get("genetic_code_cds")
    if table is None:
        return

    features = parse_annotation_gff(merged_gff)
    if not features:
        return

    sequences = {
        rec.id: rec.seq for rec in SeqIO.parse(target_fasta, "fasta")
    }
    by_gene = {}
    for feature in features:
        protein = translate_feature(sequences, feature, table)
        if not protein:
            continue
        by_gene.setdefault(feature["gene"], []).append(
            (feature["id"], protein))

    for gene, records in by_gene.items():
        lines = []
        for feature_id, protein in records:
            lines.append(f">{feature_id}")
            lines.append(protein)
        (out_dir / f"{gene}.faa").write_text("\n".join(lines) + "\n")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--merged-gff", type=Path, required=True)
    p.add_argument("--target-fasta", type=Path, required=True)
    p.add_argument("--annotation-summary", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()

    run(
        args.merged_gff, args.target_fasta, args.annotation_summary,
        args.out_dir,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
