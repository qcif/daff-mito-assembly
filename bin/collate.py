#!/usr/bin/env python3
"""COLLATE (C6) — spec §2 stage 15, §6a.5, task 42.

Reads every upstream artifact for one sample, derives a single
`sample_status` (CONSTITUTION.md principle 7 — five distinguishable
signals, never a boolean), assembles the per-sample output bundle
(`organelle_assembly.fasta`, `organelle_annotation.gff`,
`barcodes.fasta`, `metadata.json`, `report.html`), and emits
`metadata.json` — the Taxodactyl handoff record, the biosecurity audit
trail (rule 16), and the input contract for both report tiers.

`report.html` is a placeholder here — task 43_per_sample_report.md
fills it in; this task ships the data layer only.

Always exits 0 (CONSTITUTION rule 8, mirroring C2/C8's contract): a
sample whose inputs are unreadable gets a `metadata.json` recording
that fact, not a stack trace that takes its siblings down.
"""

import argparse
import json
import sys
from pathlib import Path

# Sibling modules in bin/, staged onto the container PATH by Nextflow's
# bin/ auto-staging — the same mechanism bin_target.py relies on for
# intervals/plastid_canonicalise (rule 19: reuse, don't reimplement).
from annotate_summary import build_alias_index, canonicalise

SCHEMA = "wf5/sample-metadata/v1"

# sample_status.json's gate vocabulary (task 35). Anything else fails
# closed (§3.1) rather than defaulting to "ok".
GATE_OK = "ok"
GATE_LOW_COVERAGE = "low_coverage"
GATE_FAIL = "fail"
KNOWN_GATE_STATUSES = (GATE_OK, GATE_LOW_COVERAGE, GATE_FAIL)

# Derived sample_status vocabulary (CONSTITUTION principle 7, §3.1).
STATUS_OK = "ok"
STATUS_LOW_COVERAGE = "low_coverage"
STATUS_NO_ASSEMBLY = "no_assembly"
STATUS_NO_BARCODE = "no_barcode"
STATUS_FAIL = "fail"

FULL_BUNDLE_STATUSES = (STATUS_OK, STATUS_LOW_COVERAGE, STATUS_NO_BARCODE)

BUNDLE_FULL = "full"
BUNDLE_MINIMAL = "minimal"

DIAGNOSTICS_DIR = "diagnostics"

# Reasonable defaults for tools that carry no self-reported version
# (blastn, Flye, BandageNG) — pinned per conf/containers.config, the
# single source of truth for the actual running version.
TOOL_VERSIONS = {
    "flye": "2.9.6",
    "bandage_ng": "2026.6.1",
    "blastn": "2.17.0",
}


def read_json(path):
    """Return the parsed JSON at `path`, or None if it is missing,
    empty, or malformed. Absence is data, not an error (§7)."""
    if not path:
        return None
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return None
    try:
        return json.loads(p.read_text())
    except json.JSONDecodeError as exc:
        print(
            f"WARNING: {p} is not valid JSON ({exc}) — treated as absent",
            file=sys.stderr,
        )
        return None


def read_tsv(path):
    """Return a list of dict rows from a tab-delimited file with a
    header row, or [] if the file is missing or empty."""
    if not path:
        return []
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return []
    lines = p.read_text().splitlines()
    header = lines[0].split("\t")
    return [
        dict(zip(header, line.split("\t")))
        for line in lines[1:]
        if line
    ]


def is_empty_fasta(path) -> bool:
    """A `target.fasta` that does not exist, or exists with zero bytes,
    both count as empty — BIN_TARGET's withheld-substitution case
    writes the latter (§3.2)."""
    if not path:
        return True
    p = Path(path)
    return not p.is_file() or p.stat().st_size == 0


def classify_sample(
    gate_status, contigs_selected, target_fasta_empty, n_loci_passed,
):
    """Derive the sample-level status per §3.1's dispatch matrix.

    Dispatch on a parsed, explicit allowlist; fail closed on an
    unrecognised gate status rather than defaulting to "ok" (rule 18).
    """
    if gate_status not in KNOWN_GATE_STATUSES:
        return (
            STATUS_FAIL,
            f"unrecognised coverage-gate status {gate_status!r} — "
            "failing closed to the minimal bundle",
        )

    if gate_status == GATE_FAIL:
        return STATUS_FAIL, "coverage gate soft-failed below the hard floor"

    no_assembly = (not contigs_selected) or target_fasta_empty
    if no_assembly:
        return (
            STATUS_NO_ASSEMBLY,
            "coverage gate passed but no target contig was selected "
            "(or the withheld plastid substitution left target.fasta "
            "empty — spec §3.6)",
        )

    if n_loci_passed == 0:
        return (
            STATUS_NO_BARCODE,
            "assembly and annotation exist, but no barcode locus "
            "passed validation",
        )

    if gate_status == GATE_LOW_COVERAGE:
        return (
            STATUS_LOW_COVERAGE,
            "assembled under the coverage warn floor — result is real "
            "but incomplete",
        )

    return STATUS_OK, "clean recovery"


def bundle_kind(sample_status: str) -> str:
    return (
        BUNDLE_FULL if sample_status in FULL_BUNDLE_STATUSES
        else BUNDLE_MINIMAL
    )


def kingdom_organelle(assembly_target: str) -> tuple:
    kingdom = "plant" if assembly_target.startswith("plant_") else "animal"
    organelle = assembly_target.rsplit("_", 1)[-1]
    return kingdom, organelle


def parse_assembly_info(path):
    """Flye assembly_info.txt -> list of {contig, length, coverage,
    circular} dicts, or [] if the file is missing/empty (no assembly
    was attempted)."""
    if not path or not Path(path).is_file():
        return []
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        rows.append({
            "contig": parts[0],
            "length": int(parts[1]),
            "coverage": float(parts[2]),
            "circular": parts[3].strip().upper() == "Y",
        })
    return rows


def n50(lengths) -> int:
    if not lengths:
        return None
    ordered = sorted(lengths, reverse=True)
    total = sum(ordered)
    half = total / 2
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return ordered[-1]  # pragma: no cover — unreachable, guards division


def assembly_section(assembly_info_rows, bin_metadata):
    if not assembly_info_rows and bin_metadata is None:
        return None
    lengths = [r["length"] for r in assembly_info_rows]
    return {
        "contig_count": len(assembly_info_rows) or None,
        "n50": n50(lengths),
        "total_bp": sum(lengths) if lengths else None,
        "contigs": assembly_info_rows,
        "bin_metadata": bin_metadata,
    }


def top_blast_hits(path):
    """Best hit (max bitscore) per query contig from blastn outfmt 6
    `qaccver saccver pident length qcovs evalue bitscore stitle`."""
    cols = [
        "qaccver", "saccver", "pident", "length",
        "qcovs", "evalue", "bitscore", "stitle",
    ]
    best = {}
    if not path or not Path(path).is_file():
        return []
    for line in Path(path).read_text().splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < len(cols):
            continue
        row = dict(zip(cols, fields))
        row["pident"] = float(row["pident"])
        row["length"] = int(row["length"])
        row["qcovs"] = float(row["qcovs"])
        row["evalue"] = float(row["evalue"])
        row["bitscore"] = float(row["bitscore"])
        current = best.get(row["qaccver"])
        if current is None or row["bitscore"] > current["bitscore"]:
            best[row["qaccver"]] = row
    return [best[k] for k in sorted(best)]


def homology_section(blast_tsv_path):
    hits = top_blast_hits(blast_tsv_path)
    if not hits:
        return None
    return {"top_hits": hits}


def barcodes_section(validation_tsv_path):
    rows = read_tsv(validation_tsv_path)
    if not rows:
        return None
    n_passed = sum(1 for r in rows if r.get("status") == "pass")
    return {"loci": rows, "n_passed": n_passed}


def n_loci_passed(validation_tsv_path) -> int:
    rows = read_tsv(validation_tsv_path)
    return sum(1 for r in rows if r.get("status") == "pass")


def with_canonical_names(summary: dict, gene_sets_path) -> dict:
    """Add a canonical spelling alongside every raw annotator gene name
    in `annotation_summary.json`'s cross-check buckets (§5.3 —
    `cds_crosscheck.agreed`/`miniprot_only` are flat raw-name lists;
    `annotator_only` is a list of {gene, reason} dicts). `cds_rescued`
    is left untouched: rescued features are already named by their
    canonical gene (`build_rescued_feature`), so there is no raw form
    to pair it with."""
    crosscheck = summary.get("cds_crosscheck")
    if not crosscheck or not gene_sets_path:
        return summary

    gene_sets = read_json(gene_sets_path) or {}
    gene_set = gene_sets.get("sets", {}).get(summary.get("assembly_target"))
    alias_index = build_alias_index(gene_set) if gene_set else {}

    def pair(names):
        return {
            "raw": list(names),
            "canonical": [canonicalise(n, alias_index) for n in names],
        }

    out = dict(summary)
    out_crosscheck = dict(crosscheck)
    for flat_field in ("agreed", "miniprot_only"):
        if flat_field in out_crosscheck:
            out_crosscheck[flat_field] = pair(out_crosscheck[flat_field])
    if "annotator_only" in out_crosscheck and out_crosscheck["annotator_only"]:
        out_crosscheck["annotator_only"] = [
            {**entry, "canonical": canonicalise(entry["gene"], alias_index)}
            for entry in out_crosscheck["annotator_only"]
        ]
    out["cds_crosscheck"] = out_crosscheck
    return out


def provenance_section(
    annotation_summary, genetic_code, refs_manifest, pipeline_commit,
):
    tool_versions = dict(TOOL_VERSIONS)
    if annotation_summary:
        tool_versions.update(annotation_summary.get("tool_versions") or {})
    return {
        "pipeline_commit": pipeline_commit or None,
        "reference_bundle_version": (
            (refs_manifest or {}).get("version")
        ),
        "reference_bundle_generated_at": (
            (refs_manifest or {}).get("generated_at")
        ),
        "tool_versions": tool_versions,
        "genetic_code": genetic_code,
    }


def build_metadata(args) -> dict:
    meta = read_json(args.meta_json) or {}
    status_json = read_json(args.status_json) or {}
    coverage_json = read_json(args.coverage_json) or {}
    bin_metadata = read_json(args.bin_metadata_json)
    annotation_summary = read_json(args.annotation_summary)
    if annotation_summary is not None:
        annotation_summary = with_canonical_names(
            annotation_summary, args.gene_sets)
    genetic_code = read_json(args.genetic_code_json)
    refs_manifest = read_json(args.refs_manifest)

    gate_status = status_json.get("status")
    contigs_selected = (
        (bin_metadata or {}).get("contigs_selected") or []
    )
    target_empty = is_empty_fasta(args.target_fasta)
    loci_passed = n_loci_passed(args.validation_tsv)

    sample_status, reason = classify_sample(
        gate_status, contigs_selected, target_empty, loci_passed)

    assembly_target = meta["assembly_target"]
    kingdom, organelle = kingdom_organelle(assembly_target)
    assembly_info_rows = parse_assembly_info(args.assembly_info)

    return {
        "$schema": SCHEMA,
        "sample_id": meta["sample_id"],
        "assembly_target": assembly_target,
        "kingdom": kingdom,
        "organelle": organelle,
        "sample_info": meta.get("sample_info", ""),
        "sample_type": meta.get("sample_type", ""),
        "sample_receipt_date": meta.get("sample_receipt_date", ""),
        "storage_location": meta.get("storage_location", ""),
        "sample_status": sample_status,
        "sample_status_reason": reason,
        "bundle": bundle_kind(sample_status),
        "coverage": {
            "gate": status_json or None,
            "estimate": coverage_json or None,
        },
        "assembly": assembly_section(assembly_info_rows, bin_metadata),
        "homology": homology_section(args.blast_tsv),
        "barcodes": barcodes_section(args.validation_tsv),
        "annotation": annotation_summary,
        "provenance": provenance_section(
            annotation_summary, genetic_code, refs_manifest,
            args.pipeline_commit,
        ),
    }


def copy_bytes(src, dst) -> None:
    Path(dst).write_bytes(Path(src).read_bytes())


def assemble_bundle(args, metadata: dict) -> None:
    """Write the per-sample bundle: full (assembly + annotation +
    barcodes) or minimal (metadata + report only), plus diagnostics —
    all relative to the process work directory, which `publishDir`
    copies verbatim into `outdir/<sample_id>/` (spec §0). No empty
    placeholder FASTA/GFF on the minimal path — an empty file that
    looks like a result is worse than an absent one (rule 18)."""
    diagnostics = Path(DIAGNOSTICS_DIR)

    if metadata["bundle"] == BUNDLE_FULL:
        if args.target_fasta and Path(args.target_fasta).is_file():
            copy_bytes(args.target_fasta, "organelle_assembly.fasta")
        if args.annotation_gff and Path(args.annotation_gff).is_file():
            copy_bytes(args.annotation_gff, "organelle_annotation.gff")
        if args.barcodes_fasta and Path(args.barcodes_fasta).is_file():
            copy_bytes(args.barcodes_fasta, "barcodes.fasta")

    diagnostic_files = {
        "validation.tsv": args.validation_tsv,
        "secondaries.tsv": args.secondaries_tsv,
        "bin_metadata.json": args.bin_metadata_json,
        "annotation_summary.json": args.annotation_summary,
        "organelle_map.svg": args.organelle_map_svg,
        "graph.png": args.graph_png,
    }
    has_diagnostics = any(
        p and Path(p).is_file() for p in diagnostic_files.values()
    ) or (args.plastid_isoforms and Path(args.plastid_isoforms).is_dir())
    if has_diagnostics:
        diagnostics.mkdir(parents=True, exist_ok=True)
        for name, path in diagnostic_files.items():
            if path and Path(path).is_file():
                copy_bytes(path, diagnostics / name)
        if args.plastid_isoforms and Path(args.plastid_isoforms).is_dir():
            iso_dir = diagnostics / "plastid_isoforms"
            iso_dir.mkdir(exist_ok=True)
            for child in Path(args.plastid_isoforms).iterdir():
                copy_bytes(child, iso_dir / child.name)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--meta-json", type=Path, required=True,
        help="JSON file carrying sample_id, assembly_target, and the "
             "submitter-supplied optional samplesheet columns — written "
             "by the calling Nextflow process rather than passed as "
             "individual free-text CLI flags, since sample_info etc. "
             "are unconstrained submitter text and shell-quoting them "
             "directly would be a command-injection risk")
    p.add_argument("--pipeline-commit", default="")

    p.add_argument("--status-json", type=Path, required=True)
    p.add_argument("--coverage-json", type=Path, required=True)
    p.add_argument("--bin-metadata-json", type=Path, default=None)
    p.add_argument("--assembly-info", type=Path, default=None)
    p.add_argument("--target-fasta", type=Path, default=None)
    p.add_argument("--blast-tsv", type=Path, default=None)
    p.add_argument("--barcodes-fasta", type=Path, default=None)
    p.add_argument("--secondaries-tsv", type=Path, default=None)
    p.add_argument("--validation-tsv", type=Path, default=None)
    p.add_argument("--annotation-gff", type=Path, default=None)
    p.add_argument("--annotation-summary", type=Path, default=None)
    p.add_argument("--genetic-code-json", type=Path, default=None)
    p.add_argument("--organelle-map-svg", type=Path, default=None)
    p.add_argument("--graph-png", type=Path, default=None)
    p.add_argument("--plastid-isoforms", type=Path, default=None)

    p.add_argument("--gene-sets", type=Path, default=None)
    p.add_argument("--refs-manifest", type=Path, default=None)
    p.add_argument("--schema", type=Path, default=None,
                   help="assets/sample_metadata.schema.json — validated "
                        "against if given, warned to stderr on mismatch")

    p.add_argument("--out-metadata", type=Path, default=Path("metadata.json"))
    p.add_argument("--out-report", type=Path, default=Path("report.html"))
    args = p.parse_args()

    metadata = build_metadata(args)
    assemble_bundle(args, metadata)

    if args.schema is not None and args.schema.is_file():
        try:
            import jsonschema
            jsonschema.validate(
                metadata, json.loads(args.schema.read_text()))
        except Exception as exc:  # noqa: BLE001 — never fail the sample
            print(
                f"WARNING: metadata.json failed schema validation: {exc}",
                file=sys.stderr,
            )

    args.out_metadata.write_text(json.dumps(metadata, indent=2) + "\n")
    args.out_report.touch()  # task 43_per_sample_report.md fills this in
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
