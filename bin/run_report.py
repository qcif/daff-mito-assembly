#!/usr/bin/env python3
"""RUN_REPORT (C7) — spec §2 stage 16, §6a.4, task 45.

Cross-sample join point, the last stage of the pipeline. Reads every
per-sample `metadata.json` COLLATE (C6) produced and reconciles the
sample_ids actually seen against the samplesheet: COLLATE always emits
a `metadata.json` for a sample it reaches, even on the `fail` /
`no_assembly` branches (see modules/local/collate.nf), so a
samplesheet sample_id with no matching bundle crashed somewhere
upstream of COLLATE. That is the one status value C7 derives itself —
`error` — added to the five-value `sample_status` vocabulary of
CONSTITUTION.md principle 7.

Emits `run_manifest.json` (samplesheet snapshot + hash, the full
reference-bundle manifest, pipeline commit, invocation timestamp — the
run-level half of rule 16's audit trail) and `run-report.html`,
rendered via the same Jinja machinery task 43a built for the per-sample
report (`report/run_report.py`).

Always exits 0 (CONSTITUTION rule 8/18, mirroring C6): a classification
or rendering defect here must never take the whole run's outputs down
with it — see `render_run_report` below.
"""

import argparse
import csv
import hashlib
import html
import json
import sys
import traceback
from pathlib import Path

from report.run_report import render

# Five-value vocabulary of CONSTITUTION.md principle 7, plus `error` —
# the value C7 derives itself (spec §6a.4). Order here is the display
# order for the status-breakdown table.
STATUS_ORDER = [
    "ok", "low_coverage", "no_assembly", "no_barcode", "fail", "error",
]


def read_json(path):
    """Return the parsed JSON at `path`, or None if it is missing,
    empty, or malformed. Absence is data, not an error — mirrors
    collate.py's helper of the same name."""
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


def read_samplesheet(path) -> list:
    """Return the samplesheet's rows as an ordered list of dicts — the
    reconciliation source for `error` detection, and the `rows` snapshot
    embedded in run_manifest.json."""
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh))


def sha256_of(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_bundles(metadata_paths) -> dict:
    """Every collected metadata.json, keyed by sample_id. A path that
    fails to parse, or a bundle with no sample_id, is silently dropped
    rather than raising — the sample it belonged to then reconciles as
    `error` below regardless of *why* its metadata.json didn't make it
    in, which is the honest state: C7 cannot tell a crash apart from a
    corrupt bundle, and both mean "no trustworthy record survived"."""
    bundles = {}
    for path in metadata_paths:
        metadata = read_json(path)
        sample_id = (metadata or {}).get("sample_id")
        if sample_id:
            bundles[sample_id] = metadata
    return bundles


def sample_row(sample_id: str, metadata) -> dict:
    """One row for both run_manifest.json's `samples` list and the
    run-report.html summary table (spec §6a.4's per-sample columns:
    kingdom, gate status, assembly outcome, coverage, top BLAST hit,
    panel recovery count)."""
    if metadata is None:
        return {
            "sample_id": sample_id,
            "kingdom": None,
            "assembly_target": None,
            "sample_status": "error",
            "sample_status_reason": (
                "no metadata.json was produced for this sample — the "
                "pipeline did not reach COLLATE"
            ),
            "gate_status": None,
            "coverage": None,
            "top_blast_hit": None,
            "n_barcodes_passed": None,
            "n_barcodes_total": None,
        }

    gate = (metadata.get("coverage") or {}).get("gate") or {}
    barcodes = metadata.get("barcodes") or {}
    loci = barcodes.get("loci") or []
    hits = (metadata.get("homology") or {}).get("top_hits") or []
    top_hit = max(hits, key=lambda h: h.get("bitscore", 0)) if hits else None

    return {
        "sample_id": sample_id,
        "kingdom": metadata.get("kingdom"),
        "assembly_target": metadata.get("assembly_target"),
        "sample_status": metadata.get("sample_status"),
        "sample_status_reason": metadata.get("sample_status_reason"),
        "gate_status": gate.get("status"),
        "coverage": gate.get("estimated_cov"),
        "top_blast_hit": top_hit,
        "n_barcodes_passed": barcodes.get("n_passed"),
        "n_barcodes_total": len(loci) if loci else None,
    }


def build_manifest(args) -> dict:
    bundles = load_bundles(args.metadata_json)
    sheet_rows = read_samplesheet(args.samplesheet)
    sample_ids = [
        row["sample_id"] for row in sheet_rows if row.get("sample_id")
    ]

    samples = [
        sample_row(sample_id, bundles.get(sample_id))
        for sample_id in sample_ids
    ]

    status_breakdown = {status: 0 for status in STATUS_ORDER}
    for row in samples:
        status = row["sample_status"]
        if status in status_breakdown:
            status_breakdown[status] += 1

    return {
        "$schema": "wf5/run-manifest/v1",
        "pipeline_commit": args.pipeline_commit or None,
        "invocation_timestamp": args.workflow_start,
        "samplesheet": {
            "path": str(args.samplesheet),
            "sha256": sha256_of(args.samplesheet),
            "rows": sheet_rows,
        },
        "reference_bundle": read_json(args.refs_manifest),
        "sample_count": len(samples),
        "status_breakdown": status_breakdown,
        "samples": samples,
    }


def render_run_report(args, manifest: dict) -> None:
    """Render run-report.html, never letting a rendering defect take
    the run's other outputs down with it (mirrors collate.py's
    render_bundle_report — CONSTITUTION rule 8). Unlike COLLATE's
    per-sample fallback, there is no next stage waiting on this file;
    the failure mode that matters is `run_manifest.json` staying valid
    even when the HTML does not."""
    if args.report_templates is None or args.report_static is None:
        args.out_report.touch()
        return
    try:
        render(
            manifest=manifest,
            template_dir=args.report_templates,
            static_dir=args.report_static,
            out_path=args.out_report,
        )
    except Exception:  # noqa: BLE001 — rendering must never fail the run
        tb = traceback.format_exc()
        print(f"WARNING: run-report rendering failed:\n{tb}", file=sys.stderr)
        args.out_report.write_text(
            "<html><body>"
            "<h1>Run report rendering failed</h1>"
            "<p>run_manifest.json was written successfully; the "
            "human-readable run-report.html failed to render. See the "
            "traceback below for diagnosis.</p>"
            f"<pre>{html.escape(tb)}</pre>"
            "</body></html>"
        )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--metadata-json", type=Path, nargs="+", default=[],
        help="Every collected per-sample metadata.json (C6 output)")
    p.add_argument("--samplesheet", type=Path, required=True)
    p.add_argument("--refs-manifest", type=Path, default=None,
                   help="refs/<version>/manifest.json (spec §4.4)")
    p.add_argument("--pipeline-commit", default="")
    p.add_argument("--workflow-start", default=None,
                   help="workflow.start (ISO 8601), the invocation "
                        "timestamp recorded in run_manifest.json")
    p.add_argument("--schema", type=Path, default=None,
                   help="assets/run_manifest.schema.json")
    p.add_argument("--report-templates", type=Path, default=None)
    p.add_argument("--report-static", type=Path, default=None)
    p.add_argument("--out-manifest", type=Path,
                   default=Path("run_manifest.json"))
    p.add_argument("--out-report", type=Path,
                   default=Path("run-report.html"))
    args = p.parse_args()

    manifest = build_manifest(args)

    if args.schema is not None and args.schema.is_file():
        try:
            import jsonschema
            jsonschema.validate(
                manifest, json.loads(args.schema.read_text()))
        except Exception as exc:  # noqa: BLE001 — never fail the run
            print(
                f"WARNING: run_manifest.json failed schema validation: "
                f"{exc}",
                file=sys.stderr,
            )

    args.out_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    render_run_report(args, manifest)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
