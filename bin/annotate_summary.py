#!/usr/bin/env python3
"""ANNOTATE (C8) — spec §2 stage 13a, §2.2 C8, task 31.

Merges ``MINIPROT_CDS``'s ``cds.gff`` with a specialist annotator's
non-CDS (tRNA/rRNA) features into one GFF3, keyed by contig name, and
writes ``annotation_summary.json``. Today the only specialist
annotator is MITOS2, run on ``animal_mt`` by the calling Nextflow
process; the plant branches have none yet (task 32, task 33) and stay
CDS-only.

**CDS features are copied verbatim from ``cds.gff``** — coordinates
are never adjusted, re-derived or re-sorted. The comprehensive protein
panel (task 29) carries several representative sequences per gene, so
more than one ``cds.gff`` hit commonly maps to the same genomic locus;
overlapping hits for the same gene are clustered and only the
highest-identity hit per cluster is emitted (mirroring
``validate_barcodes.py``'s locus clustering). Non-overlapping clusters
for the same gene are genuine distinct loci (inverted-repeat
duplication, plant-mt repeats) and are never collapsed.

**MITOS2's own protein-coding calls are never emitted as features.**
They are compared against the ``cds.gff`` winners as an independent
cross-check (``cds_crosscheck``) and discarded otherwise — emitting
both would let the report show a MITOS2 ``cox1`` at different
coordinates from the miniprot ``cox1`` the barcode was cut from,
exactly the incoherence the unified miniprot pass (task 30) removed.

Stdlib-only: this module runs inside the MITOS2 biocontainer, which
carries no ``biopython``/``pandas`` (task 31 §2).

Always exits 0 (CONSTITUTION rule 8) — a poor or missing annotation is
supplementary information, not a pipeline error.
"""

import argparse
import json
import re
import sys
from pathlib import Path

GFF_COLUMNS = [
    "seqid", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes",
]

STATUS_NO_ASSEMBLY = "no_assembly"
STATUS_OK_CDS_ONLY = "ok_cds_only"
STATUS_ANNOTATOR_FAILED = "annotator_failed"
STATUS_OK = "ok"
STATUS_NO_FEATURES = "no_features"

SCHEMA = "wf5/annotation-summary/v1"

MINIPROT_VERSION = "0.18"
MITOS_VERSION = "2.1.10"


def parse_attrs(field: str) -> dict:
    attrs = {}
    for kv in field.strip().split(";"):
        if not kv or "=" not in kv:
            continue
        key, value = kv.split("=", 1)
        attrs[key] = value
    return attrs


def parse_gff_line(line: str, line_no: int, path: Path):
    """Parse one tab-delimited GFF3 line into a dict, or None + a
    stderr warning if malformed (spec principle 8 — never abort)."""
    fields = line.split("\t")
    if len(fields) != 9:
        print(
            f"WARNING: skipping malformed GFF line {line_no} in "
            f"{path}: expected 9 columns, got {len(fields)}",
            file=sys.stderr,
        )
        return None
    try:
        seqid, source, ftype, start, end, score, strand, phase, attr = \
            fields
        record = dict(zip(GFF_COLUMNS, fields))
        record["start"] = int(start)
        record["end"] = int(end)
        record["attrs"] = parse_attrs(attr)
        record["raw"] = line
        return record
    except (ValueError, IndexError) as exc:
        print(
            f"WARNING: skipping malformed GFF line {line_no} in "
            f"{path}: {exc}",
            file=sys.stderr,
        )
        return None


# ── cds.gff (miniprot) parsing ──────────────────────────────────────

def parse_cds_gff(path: Path) -> list:
    """Parse mRNA + CDS rows into one record per miniprot hit.

    Mirrors ``validate_barcodes.parse_gff`` but keeps the raw mRNA
    line too (needed to reproduce it verbatim in the merged GFF) and
    does not need CIGAR/PAF context, since C8 does no translation.
    """
    mrna = {}
    with open(path) as fh:
        for line_no, raw_line in enumerate(fh, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            record = parse_gff_line(line, line_no, path)
            if record is None:
                continue
            if record["type"] == "mRNA":
                attrs = record["attrs"]
                if "ID" not in attrs or "Target" not in attrs:
                    print(
                        f"WARNING: skipping malformed GFF line "
                        f"{line_no} in {path}: mRNA missing ID/Target",
                        file=sys.stderr,
                    )
                    continue
                target = attrs["Target"].split()
                protein_id = target[0]
                gene = protein_id.rsplit("_", 1)[-1]
                mrna[attrs["ID"]] = {
                    "seqid": record["seqid"],
                    "strand": record["strand"],
                    "gene": gene,
                    "identity": float(attrs.get("Identity", 0.0)),
                    "mrna_line": line,
                    "cds_lines": [],
                    "cds": [],
                }
            elif record["type"] == "CDS":
                parent = record["attrs"].get("Parent")
                if parent in mrna:
                    mrna[parent]["cds_lines"].append(line)
                    mrna[parent]["cds"].append(
                        (record["start"], record["end"]))
    return [rec for rec in mrna.values() if rec["cds"]]


def cluster_cds_by_gene(records: list) -> dict:
    """Group cds.gff records by gene symbol, then merge genomically
    overlapping hits (redundant panel representatives of the same
    locus), keeping the highest-identity hit per cluster.
    Non-overlapping hits are genuine distinct loci and stay separate
    (task 31 §4 point 6 — never de-duplicate by gene name)."""
    by_gene = {}
    for rec in records:
        rec["start"] = min(s for s, _e in rec["cds"])
        rec["end"] = max(e for _s, e in rec["cds"])
        by_gene.setdefault(rec["gene"], []).append(rec)

    winners_by_gene = {}
    for gene, recs in by_gene.items():
        recs = sorted(recs, key=lambda r: (r["seqid"], r["start"]))
        clusters = []
        for rec in recs:
            if (clusters
                    and clusters[-1][-1]["seqid"] == rec["seqid"]
                    and rec["start"]
                    <= max(r["end"] for r in clusters[-1])):
                clusters[-1].append(rec)
            else:
                clusters.append([rec])
        winners_by_gene[gene] = [
            max(cluster, key=lambda r: r["identity"])
            for cluster in clusters
        ]
    return winners_by_gene


# ── MITOS2 result.gff parsing ───────────────────────────────────────

def find_mitos_result_gffs(mitos_dir: Path) -> list:
    """MITOS2 writes ``result.gff`` directly under ``--outdir`` for a
    single-record input, or under one numbered subdirectory per input
    FASTA record (``<outdir>/<i>/result.gff``) for multi-record input.
    Either way, the seqid column already carries the real input record
    id (``mitos.gfffile.gffwriter`` writes the accession it is passed,
    which is the record's own id) — no index-based remapping needed,
    just gather every ``result.gff`` found."""
    return sorted(mitos_dir.rglob("result.gff"))


def parse_mitos_gff(mitos_dir: Path) -> tuple:
    """Returns (non_cds_records, cds_call_records).

    MITOS2 emits three rows per RNA gene (``ncRNA_gene``, ``tRNA`` or
    ``rRNA``, ``exon``) and two per protein gene (``gene``, ``exon``).
    Only the summary row per gene is used — ``tRNA``/``rRNA`` rows are
    the non-CDS features emitted into the merged GFF; ``gene`` rows are
    MITOS2's own CDS calls, used only for the cross-check (§3/§4).
    """
    non_cds = []
    cds_calls = []
    for gff_path in find_mitos_result_gffs(mitos_dir):
        with open(gff_path) as fh:
            for line_no, raw_line in enumerate(fh, 1):
                line = raw_line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                record = parse_gff_line(line, line_no, gff_path)
                if record is None:
                    continue
                if record["type"] in ("tRNA", "rRNA"):
                    record["gene"] = record["attrs"].get("Name", "")
                    non_cds.append(record)
                elif record["type"] == "gene":
                    record["gene"] = record["attrs"].get("Name", "")
                    cds_calls.append(record)
    return non_cds, cds_calls


# ── gene-name normalisation + completeness ──────────────────────────

ANTICODON_RE = re.compile(r"\(.*?\)")


def normalise_gene_symbol(raw: str) -> str:
    """Fold a gene symbol for cross-source comparison: uppercase,
    strip a MITOS2 anticodon parenthetical (``trnF(gaa)`` ->
    ``TRNF``), strip an ``MT-`` locus prefix, drop spaces/dashes/
    underscores. The raw name in GFF attributes is left untouched —
    only comparison uses this folded form (task 31 §4 point 4)."""
    s = ANTICODON_RE.sub("", raw).upper()
    if s.startswith("MT-"):
        s = s[3:]
    return re.sub(r"[\s_-]", "", s)


def build_alias_index(gene_set: dict) -> dict:
    """normalised alias -> canonical protein_coding gene name."""
    index = {}
    aliases = gene_set.get("protein_coding_aliases", {})
    for canonical in gene_set["protein_coding"]:
        index[normalise_gene_symbol(canonical)] = canonical
        for alias in aliases.get(canonical, []):
            index[normalise_gene_symbol(alias)] = canonical
    return index


def protein_coding_completeness(winners_by_gene: dict, gene_set: dict):
    alias_index = build_alias_index(gene_set)
    found = set()
    for gene in winners_by_gene:
        canonical = alias_index.get(normalise_gene_symbol(gene))
        if canonical is not None:
            found.add(canonical)
    canonical_set = set(gene_set["protein_coding"])
    missing = sorted(canonical_set - found)
    fraction = (
        len(found) / len(canonical_set) if canonical_set else None
    )
    return fraction, sorted(found), missing


# ── cross-check ──────────────────────────────────────────────────────

def overlaps(a: dict, b: dict) -> bool:
    return (
        a["seqid"] == b["seqid"]
        and a["start"] <= b["end"]
        and b["start"] <= a["end"]
    )


def canonicalise(raw: str, alias_index: dict) -> str:
    """Fold `raw` to its gene-sets canonical symbol where a mapping
    exists (so e.g. MITOS2's ``cob`` and miniprot's ``CYTB`` compare
    equal), else fall back to the plain normalised form."""
    normalised = normalise_gene_symbol(raw)
    return alias_index.get(normalised, normalised)


def cds_crosscheck(
    winners_by_gene: dict, mitos_cds_calls: list, alias_index: dict,
) -> dict:
    miniprot_by_norm = {}
    for gene, winners in winners_by_gene.items():
        miniprot_by_norm.setdefault(
            canonicalise(gene, alias_index), []).extend(winners)

    annotator_by_norm = {}
    for call in mitos_cds_calls:
        annotator_by_norm.setdefault(
            canonicalise(call["gene"], alias_index), []).append(call)

    agreed, annotator_only, coordinate_conflicts = [], [], []
    for norm_name, calls in annotator_by_norm.items():
        miniprot_hits = miniprot_by_norm.get(norm_name)
        if not miniprot_hits:
            annotator_only.append(calls[0]["gene"])
            continue
        if any(
            overlaps(call, hit)
            for call in calls for hit in miniprot_hits
        ):
            agreed.append(calls[0]["gene"])
        else:
            coordinate_conflicts.append(calls[0]["gene"])

    miniprot_only = [
        winners[0]["gene"]
        for norm_name, winners in miniprot_by_norm.items()
        if norm_name not in annotator_by_norm
    ]

    return {
        "agreed": sorted(agreed),
        "miniprot_only": sorted(miniprot_only),
        "annotator_only": sorted(annotator_only),
        "coordinate_conflicts": sorted(coordinate_conflicts),
    }


# ── output writers ───────────────────────────────────────────────────

def write_gff(out_path: Path, winners_by_gene: dict, non_cds: list):
    lines = ["##gff-version 3"]
    cds_blocks = [
        winner
        for winners in winners_by_gene.values()
        for winner in winners
    ]
    cds_blocks.sort(key=lambda w: (w["seqid"], w["start"]))
    for winner in cds_blocks:
        lines.append(winner["mrna_line"])
        lines.extend(winner["cds_lines"])

    non_cds_sorted = sorted(
        non_cds, key=lambda r: (r["seqid"], r["start"]))
    for record in non_cds_sorted:
        lines.append(record["raw"])

    out_path.write_text("\n".join(lines) + "\n")


def run(
    cds_gff: Path, sample_id: str, assembly_target: str,
    gene_sets_path: Path, mitos_dir, genetic_code_annotate,
    genetic_code_barcodes, reference_data, annotator_exit, out_gff: Path,
    out_summary: Path,
) -> None:
    gene_sets = json.loads(gene_sets_path.read_text())["sets"]
    gene_set = gene_sets.get(assembly_target)

    cds_bytes = cds_gff.stat().st_size if cds_gff.exists() else 0

    if cds_bytes == 0:
        write_gff(out_gff, {}, [])
        summary = {
            "$schema": SCHEMA,
            "sample_id": sample_id,
            "assembly_target": assembly_target,
            "status": STATUS_NO_ASSEMBLY,
            "reason": "empty cds.gff — no target contig assembled upstream",
            "cds_source": "miniprot",
            "non_cds_source": None,
            "annotator_exit_code": None,
            "tool_versions": {"miniprot": MINIPROT_VERSION},
            "reference_data": None,
            "genetic_code_annotate": None,
            "genetic_code_barcodes": None,
            "genetic_code_agreement": None,
            "contigs_annotated": [],
            "feature_counts": {
                "gene": 0, "CDS": 0, "tRNA": 0, "rRNA": 0,
            },
            "protein_coding_completeness": None,
            "protein_coding_genes_missing": None,
            "cds_crosscheck": None,
            "canonical_gene_set": gene_set["name"] if gene_set else None,
        }
        out_summary.write_text(json.dumps(summary, indent=2) + "\n")
        return

    records = parse_cds_gff(cds_gff)
    winners_by_gene = cluster_cds_by_gene(records)
    n_cds = sum(len(w) for w in winners_by_gene.values())

    non_cds = []
    mitos_cds_calls = []
    non_cds_source = None
    found_result_gffs = []
    tool_versions = {"miniprot": MINIPROT_VERSION}
    if mitos_dir is not None:
        non_cds_source = "mitos2"
        tool_versions["mitos2"] = MITOS_VERSION
        non_cds, mitos_cds_calls = parse_mitos_gff(mitos_dir)
        found_result_gffs = find_mitos_result_gffs(mitos_dir)

    n_trna = sum(1 for r in non_cds if r["type"] == "tRNA")
    n_rrna = sum(1 for r in non_cds if r["type"] == "rRNA")

    if n_cds == 0 and not non_cds:
        status = STATUS_NO_FEATURES
        reason = "no features parsed from cds.gff or the non-CDS annotator"
    elif non_cds_source is None:
        status = STATUS_OK_CDS_ONLY
        reason = (
            "no non-CDS annotator for this assembly_target — "
            "see task 32 (plant_pt) / task 33 (plant_mt)"
        )
    elif not non_cds:
        status = STATUS_ANNOTATOR_FAILED
        exit_text = (
            f"exited {annotator_exit}" if annotator_exit is not None
            else "exit code unknown"
        )
        outcome = (
            "wrote no result.gff" if not found_result_gffs
            else "wrote result.gff with no tRNA/rRNA features"
        )
        reason = (
            f"{non_cds_source} {exit_text} and {outcome} "
            f"under {mitos_dir.name}/"
        )
    else:
        status = STATUS_OK
        reason = None

    completeness, found, missing = (
        protein_coding_completeness(winners_by_gene, gene_set)
        if gene_set else (None, None, None)
    )

    crosscheck = (
        cds_crosscheck(
            winners_by_gene, mitos_cds_calls,
            build_alias_index(gene_set) if gene_set else {})
        if non_cds_source is not None else None
    )

    contigs = sorted({
        w["seqid"] for winners in winners_by_gene.values()
        for w in winners
    } | {r["seqid"] for r in non_cds})

    genetic_code_agreement = (
        genetic_code_annotate == genetic_code_barcodes
        if genetic_code_annotate is not None
        and genetic_code_barcodes is not None
        else None
    )

    write_gff(out_gff, winners_by_gene, non_cds)

    summary = {
        "$schema": SCHEMA,
        "sample_id": sample_id,
        "assembly_target": assembly_target,
        "status": status,
        "reason": reason,
        "cds_source": "miniprot",
        "non_cds_source": non_cds_source,
        "annotator_exit_code": annotator_exit,
        "tool_versions": tool_versions,
        "reference_data": reference_data,
        "genetic_code_annotate": genetic_code_annotate,
        "genetic_code_barcodes": genetic_code_barcodes,
        "genetic_code_agreement": genetic_code_agreement,
        "contigs_annotated": contigs,
        "feature_counts": {
            "gene": n_cds + n_trna + n_rrna,
            "CDS": n_cds,
            "tRNA": n_trna,
            "rRNA": n_rrna,
        },
        "protein_coding_completeness": completeness,
        "protein_coding_genes_missing": missing,
        "cds_crosscheck": crosscheck,
        "canonical_gene_set": gene_set["name"] if gene_set else None,
    }
    out_summary.write_text(json.dumps(summary, indent=2) + "\n")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cds-gff", type=Path, required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--assembly-target", required=True)
    p.add_argument("--gene-sets", type=Path, required=True)
    p.add_argument("--mitos-dir", type=Path, default=None)
    p.add_argument("--genetic-code-annotate", default="",
                   help="NCBI genetic-code table MITOS2 was run with")
    p.add_argument(
        "--genetic-code-barcodes", default="",
        help="NCBI genetic-code table EXTRACT_BARCODES (C5) is "
             "configured to try for this assembly_target — a config "
             "proxy for C5's actual per-run clade-trial choice, which "
             "ANNOTATE does not consume (task 31 §3)")
    p.add_argument("--reference-data", default="",
                   help="Non-CDS annotator reference-data tag, "
                        "e.g. refseq89m")
    p.add_argument(
        "--annotator-exit", type=int, default=None,
        help="Non-CDS annotator's process exit code, e.g. MITOS2's "
             "$MITOS_EXIT — omitted when no annotator was configured")
    p.add_argument("--out-gff", type=Path, required=True)
    p.add_argument("--out-summary", type=Path, required=True)
    args = p.parse_args()

    genetic_code_annotate = (
        int(args.genetic_code_annotate)
        if args.genetic_code_annotate else None
    )
    genetic_code_barcodes = (
        int(args.genetic_code_barcodes)
        if args.genetic_code_barcodes else None
    )
    reference_data = args.reference_data or None

    run(
        args.cds_gff, args.sample_id, args.assembly_target,
        args.gene_sets, args.mitos_dir, genetic_code_annotate,
        genetic_code_barcodes, reference_data, args.annotator_exit,
        args.out_gff, args.out_summary,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
