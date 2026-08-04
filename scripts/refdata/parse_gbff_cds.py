#!/usr/bin/env python3
"""
Extract CDS translations from a RefSeq organelle GenBank flat file.

Streams a (concatenated, uncompressed) RefSeq plastid or mitochondrion
`.gbff` file record by record and emits one JSON line per CDS feature
that carries a `/gene` and a `/translation` qualifier:

    {"accession": "NC_012920.1", "organism": "Homo sapiens",
     "gene": "COX1", "translation": "MFAD..."}

Runs inside the neoformit/daff-wf5-scripts:test image (has biopython
pinned) rather than on the host — see scripts/pytest.sh for the same
pattern. Requires only biopython; no other third-party deps.

Usage (inside the container):
    python3 parse_gbff_cds.py \
        --gbff       refseq_pt.gbff \
        --out        pt_cds.jsonl \
        [--taxonomy-filter Viridiplantae]
"""

import argparse
import json
import sys

from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gbff", required=True)
    p.add_argument("--out", required=True)
    p.add_argument(
        "--taxonomy-filter", default="",
        help=(
            "optional lineage term (e.g. 'Metazoa', 'Viridiplantae') — "
            "only records whose GenBank taxonomy lineage contains this "
            "term are emitted. Used to kingdom-split the mitochondrion "
            "GBFF, which (unlike plastid) is not pre-split by kingdom."
        ),
    )
    return p.parse_args()


def main():
    args = parse_args()

    n_records = 0
    n_cds = 0
    with open(args.out, "w") as out_fh:
        for record in SeqIO.parse(args.gbff, "genbank"):
            n_records += 1
            if args.taxonomy_filter:
                lineage = record.annotations.get("taxonomy", [])
                if args.taxonomy_filter not in lineage:
                    continue
            organism = record.annotations.get("organism", "")
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                genes = feature.qualifiers.get("gene", [])
                translations = feature.qualifiers.get("translation", [])
                if not genes or not translations:
                    continue
                out_fh.write(json.dumps({
                    "accession": record.id,
                    "organism": organism,
                    "gene": genes[0],
                    "translation": translations[0],
                }) + "\n")
                n_cds += 1
            if n_records % 500 == 0:
                print(
                    f"[parse_gbff_cds] {n_records} records, "
                    f"{n_cds} CDS so far",
                    file=sys.stderr,
                )

    print(
        f"[parse_gbff_cds] done: {n_records} records, {n_cds} CDS -> "
        f"{args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
