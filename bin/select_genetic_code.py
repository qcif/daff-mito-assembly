#!/usr/bin/env python3
"""SELECT_GENETIC_CODE — spec §2 stage 12, §2.2, task 38 §2/§3.

`animal_mt` has no single correct NCBI translation table (vertebrate
table 2 vs. invertebrate table 5 differ on AGA/AGG/TGA), so
`MINIPROT_CDS` runs one miniprot pass per configured table and this
component picks between the resulting candidate ``cds.gff`` files.

**Selection criterion**: sum the best per-gene alignment score
(``mRNA`` column 6, miniprot's ``AS`` chaining score) across the
panel; ties broken on the number of distinct genes recovered, then on
trial order (the order ``--candidate`` was given on the command
line). Frameshift/StopCodon counts are deliberately not used — they
are symptoms of the wrong table, not an independent signal, and using
them here would couple the selector to miniprot's flag semantics
(task 38 §2).

Only the winning candidate's GFF is kept; the losing candidate(s) are
discarded. ``genetic_code.json`` records every candidate's score and
gene count as the audit trail for a decision made under uncertainty
(CONSTITUTION rule 16).

plant_pt/plant_mt (single-table targets) never call this component —
`MINIPROT_CDS` passes their one configured ``-T`` directly and skips
selection entirely (task 38 §3).
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

SCHEMA = "wf5/genetic-code-selection/v1"

CRITERION = (
    "sum of best per-gene miniprot alignment score (mRNA column 6) "
    "across the panel; ties broken on distinct gene count, then on "
    "the order candidates were given"
)


def gene_from_target(target: str) -> str:
    """``NC_000000_COX1 1 100`` -> ``COX1`` (mirrors
    ``annotate_summary.parse_cds_gff``'s gene derivation)."""
    protein_id = target.split()[0]
    return protein_id.rsplit("_", 1)[-1]


def parse_attrs(field: str) -> dict:
    return dict(
        kv.split("=", 1) for kv in field.strip().split(";") if "=" in kv
    )


def parse_best_scores(path: Path) -> dict:
    """Return {gene: best_alignment_score} from a miniprot GFF3's
    mRNA rows. Malformed or unusable rows are skipped with a stderr
    warning rather than aborting (CONSTITUTION rule 8)."""
    best = {}
    with open(path) as fh:
        for line_no, raw_line in enumerate(fh, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9 or fields[2] != "mRNA":
                continue
            try:
                score = float(fields[5])
            except ValueError:
                print(
                    f"WARNING: skipping malformed mRNA line {line_no} "
                    f"in {path}: non-numeric score {fields[5]!r}",
                    file=sys.stderr,
                )
                continue
            attrs = parse_attrs(fields[8])
            target = attrs.get("Target")
            if not target:
                print(
                    f"WARNING: skipping mRNA line {line_no} in {path}: "
                    "no Target attribute",
                    file=sys.stderr,
                )
                continue
            gene = gene_from_target(target)
            if score > best.get(gene, float("-inf")):
                best[gene] = score
    return best


def select(candidates: list) -> tuple:
    """candidates: [(table, path), ...] in trial order.

    Returns (winner_record, all_records). `winner_record` is the
    element of `all_records` chosen by the criterion above; ties
    resolve to the earliest-listed candidate (`max` keeps the first
    of equal elements)."""
    records = []
    for table, path in candidates:
        best_scores = parse_best_scores(path)
        records.append({
            "table": table,
            "gff": str(path),
            "total_score": sum(best_scores.values()),
            "n_genes": len(best_scores),
        })
    winner = max(records, key=lambda r: (r["total_score"], r["n_genes"]))
    return winner, records


def parse_candidate_arg(raw: str) -> tuple:
    table_str, path_str = raw.split(":", 1)
    return int(table_str), Path(path_str)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--candidate", action="append", required=True, metavar="TABLE:PATH",
        help="NCBI genetic-code table and its candidate cds.gff, e.g. "
             "2:cds.2.gff — repeat in trial order")
    p.add_argument("--out-gff", type=Path, required=True)
    p.add_argument("--out-json", type=Path, required=True)
    args = p.parse_args()

    candidates = [parse_candidate_arg(c) for c in args.candidate]
    winner, records = select(candidates)

    shutil.copyfile(winner["gff"], args.out_gff)

    result = {
        "$schema": SCHEMA,
        "selected_table": winner["table"],
        "criterion": CRITERION,
        "candidates": records,
    }
    args.out_json.write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
