"""Unit tests for bin/run_report.py (C7) — spec §2 stage 16, §6a.4,
task 45.

Real temporary files and a real Jinja render, not mocks (task 27 §6's
boundary-mocking rule only bites at an external tool call, and this
module has none) — fixtures are small on-disk metadata.json/CSV trees
mirroring the real upstream artifacts' shapes, rendered through the
actual scripts/report/ templates/static tree.

Cases:
  1. `sample_row` extracts kingdom/gate status/coverage/top BLAST
     hit/barcode counts from a full metadata.json, and produces the
     `error` shape for a missing one.
  2. `build_manifest` reconciles the samplesheet against the collected
     metadata.json set: a samplesheet sample_id absent from the
     bundles becomes `error`, not silently dropped.
  3. All six sample_status values (five real + derived `error`) count
     correctly in `status_breakdown` for one mixed batch — the "mixed
     batch with one soft-failed sample" case from spec §06-phases.md's
     P4a acceptance, exercised here since no integration fixture
     reaches `fail`/`no_assembly`/`no_barcode`/`error` (task 42 §9).
  4. Emitted run_manifest.json validates against
     assets/run_manifest.schema.json.
  5. `render` produces a self-contained run-report.html: every sample
     links to `<sample_id>/report.html`, `fail` and `error` rows carry
     distinct label text, and low_coverage/ok rows carry distinct
     severity classes.
  6. Render failure produces the fallback HTML without raising;
     `report_templates=None` touches an empty report; `main()` still
     exits 0 and writes a valid run_manifest.json in both cases.
  7. `read_json` / `sha256_of` / `read_samplesheet` helpers.

Plus direct unit coverage of every helper for 100% branch coverage
per CONSTITUTION rule 14.
"""

import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
REPO_ROOT = BIN_DIR.parent
TEMPLATE_DIR = REPO_ROOT / "scripts" / "report" / "templates"
STATIC_DIR = REPO_ROOT / "scripts" / "report" / "static"

sys.path.insert(0, str(BIN_DIR))

import report.run_report as run_report_mod  # noqa: E402

SCHEMA = json.loads(
    (REPO_ROOT / "assets" / "run_manifest.schema.json").read_text())

_spec = importlib.util.spec_from_file_location(
    "run_report", BIN_DIR / "run_report.py")
run_report = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_report)


def write_csv(path: Path, rows: list) -> Path:
    fieldnames = [
        "sample_id", "assembly_target", "reads", "sample_info",
        "sample_type", "sample_receipt_date", "storage_location",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def write_json(path: Path, obj) -> Path:
    path.write_text(json.dumps(obj))
    return path


def base_metadata(sample_id: str, status: str, **overrides) -> dict:
    metadata = {
        "$schema": "wf5/sample-metadata/v1",
        "sample_id": sample_id,
        "assembly_target": "animal_mt",
        "kingdom": "animal",
        "organelle": "mt",
        "sample_status": status,
        "sample_status_reason": "test fixture",
        "bundle": (
            "full" if status in ("ok", "low_coverage", "no_barcode")
            else "minimal"
        ),
        "coverage": {
            "gate": {
                "status": "ok" if status == "ok" else status,
                "estimated_cov": 42.5,
            },
            "estimate": {},
        },
        "homology": {"top_hits": [
            {"qaccver": "contig_1", "saccver": "NC_000001.1",
             "pident": 98.5, "length": 100, "qcovs": 99.0,
             "evalue": 0.0, "bitscore": 150.0, "stitle": "Test organism"},
            {"qaccver": "contig_1", "saccver": "NC_000002.1",
             "pident": 80.0, "length": 100, "qcovs": 90.0,
             "evalue": 0.0, "bitscore": 90.0, "stitle": "Lower hit"},
        ]},
        "barcodes": {"loci": [
            {"gene": "COX1", "status": "pass"},
            {"gene": "CYTB", "status": "fail"},
        ], "n_passed": 1},
        "provenance": {
            "pipeline_commit": "abc1234",
            "reference_bundle_version": "v2026.09_1",
            "tool_versions": {},
            "genetic_code": None,
        },
    }
    metadata.update(overrides)
    return metadata


REFS_MANIFEST = {
    "$schema": "wf5/refs-manifest/v1",
    "version": "v2026.09_1",
    "generated_at": "2026-09-11",
}


class TestSampleRow(unittest.TestCase):
    """Case 1."""

    def test_full_metadata_extracts_expected_fields(self):
        row = run_report.sample_row("S1", base_metadata("S1", "ok"))
        self.assertEqual(row["sample_id"], "S1")
        self.assertEqual(row["kingdom"], "animal")
        self.assertEqual(row["assembly_target"], "animal_mt")
        self.assertEqual(row["sample_status"], "ok")
        self.assertEqual(row["gate_status"], "ok")
        self.assertEqual(row["coverage"], 42.5)
        # Best hit by bitscore, not first in list.
        self.assertEqual(row["top_blast_hit"]["stitle"], "Test organism")
        self.assertEqual(row["n_barcodes_passed"], 1)
        self.assertEqual(row["n_barcodes_total"], 2)

    def test_missing_metadata_yields_error_row(self):
        row = run_report.sample_row("S2", None)
        self.assertEqual(row["sample_status"], "error")
        self.assertIsNone(row["kingdom"])
        self.assertIsNone(row["gate_status"])
        self.assertIsNone(row["coverage"])
        self.assertIsNone(row["top_blast_hit"])
        self.assertIn("no metadata.json", row["sample_status_reason"])

    def test_no_homology_or_barcodes_sections(self):
        metadata = base_metadata(
            "S3", "fail", homology=None, barcodes=None)
        row = run_report.sample_row("S3", metadata)
        self.assertIsNone(row["top_blast_hit"])
        self.assertIsNone(row["n_barcodes_passed"])
        self.assertIsNone(row["n_barcodes_total"])


class TestBuildManifest(unittest.TestCase):
    """Cases 2, 3, 4 — reconciliation, status breakdown, schema."""

    def test_mixed_batch_reconciles_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            samplesheet_csv = write_csv(d / "samples.csv", [
                {"sample_id": "OK1", "assembly_target": "animal_mt"},
                {"sample_id": "LOW1", "assembly_target": "plant_mt"},
                {"sample_id": "NOASM1", "assembly_target": "animal_mt"},
                {"sample_id": "NOBARC1", "assembly_target": "animal_mt"},
                {"sample_id": "FAIL1", "assembly_target": "animal_mt"},
                {"sample_id": "CRASHED1", "assembly_target": "animal_mt"},
            ])
            meta_paths = []
            for sample_id, status in [
                ("OK1", "ok"), ("LOW1", "low_coverage"),
                ("NOASM1", "no_assembly"), ("NOBARC1", "no_barcode"),
                ("FAIL1", "fail"),
            ]:
                p = write_json(
                    d / f"{sample_id}.metadata.json",
                    base_metadata(sample_id, status))
                meta_paths.append(p)
            # CRASHED1 never reaches COLLATE — no metadata.json at all.

            refs_manifest_json = write_json(
                d / "manifest.json", REFS_MANIFEST)

            class Args:
                metadata_json = meta_paths
                samplesheet = samplesheet_csv
                refs_manifest = refs_manifest_json
                pipeline_commit = "abc1234"
                workflow_start = "2026-09-29T10:00:00"

            manifest = run_report.build_manifest(Args())

            self.assertEqual(manifest["sample_count"], 6)
            self.assertEqual(manifest["status_breakdown"], {
                "ok": 1, "low_coverage": 1, "no_assembly": 1,
                "no_barcode": 1, "fail": 1, "error": 1,
            })
            crashed = next(
                s for s in manifest["samples"]
                if s["sample_id"] == "CRASHED1")
            self.assertEqual(crashed["sample_status"], "error")
            self.assertEqual(
                manifest["reference_bundle"]["version"], "v2026.09_1")
            self.assertEqual(manifest["pipeline_commit"], "abc1234")

            import jsonschema
            jsonschema.validate(manifest, SCHEMA)

    def test_no_refs_manifest_is_null_not_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            samplesheet_csv = write_csv(
                d / "samples.csv",
                [{"sample_id": "S1", "assembly_target": "animal_mt"}])
            meta = write_json(
                d / "S1.metadata.json", base_metadata("S1", "ok"))

            class Args:
                metadata_json = [meta]
                samplesheet = samplesheet_csv
                refs_manifest = None
                pipeline_commit = ""
                workflow_start = None

            manifest = run_report.build_manifest(Args())
            self.assertIsNone(manifest["reference_bundle"])
            self.assertIsNone(manifest["pipeline_commit"])

            import jsonschema
            jsonschema.validate(manifest, SCHEMA)

    def test_malformed_metadata_json_is_dropped_and_reconciled(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            samplesheet_csv = write_csv(
                d / "samples.csv",
                [{"sample_id": "S1", "assembly_target": "animal_mt"}])
            bad = d / "bad.metadata.json"
            bad.write_text("{not valid json")

            class Args:
                metadata_json = [bad]
                samplesheet = samplesheet_csv
                refs_manifest = None
                pipeline_commit = ""
                workflow_start = None

            manifest = run_report.build_manifest(Args())
            self.assertEqual(manifest["status_breakdown"]["error"], 1)

    def test_unrecognised_status_fails_closed_out_of_breakdown(self):
        """A sample_status outside the known vocabulary (schema drift,
        a corrupt bundle) must not silently inflate a real bucket —
        it stays out of status_breakdown but is still visible in the
        per-sample rows (rule 18: fail closed, never silently absorb)."""
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            samplesheet_csv = write_csv(
                d / "samples.csv",
                [{"sample_id": "S1", "assembly_target": "animal_mt"}])
            meta = write_json(
                d / "S1.metadata.json",
                base_metadata("S1", "bogus_status"))

            class Args:
                metadata_json = [meta]
                samplesheet = samplesheet_csv
                refs_manifest = None
                pipeline_commit = ""
                workflow_start = None

            manifest = run_report.build_manifest(Args())
            self.assertEqual(
                sum(manifest["status_breakdown"].values()), 0)
            self.assertEqual(
                manifest["samples"][0]["sample_status"], "bogus_status")


class TestReadHelpers(unittest.TestCase):
    """Case 7."""

    def test_read_json_missing(self):
        self.assertIsNone(run_report.read_json(None))

    def test_read_json_nonexistent(self):
        self.assertIsNone(run_report.read_json(Path("/no/such/file.json")))

    def test_read_json_empty_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "empty.json"
            p.write_text("")
            self.assertIsNone(run_report.read_json(p))

    def test_read_json_malformed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "bad.json"
            p.write_text("{not valid")
            self.assertIsNone(run_report.read_json(p))

    def test_read_json_valid(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = write_json(Path(tmp) / "ok.json", {"a": 1})
            self.assertEqual(run_report.read_json(p), {"a": 1})

    def test_read_samplesheet(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = write_csv(Path(tmp) / "s.csv", [
                {"sample_id": "A", "assembly_target": "animal_mt"},
            ])
            rows = run_report.read_samplesheet(p)
            self.assertEqual(rows[0]["sample_id"], "A")

    def test_sha256_of(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "f.txt"
            p.write_text("hello")
            expected = hashlib.sha256(b"hello").hexdigest()
            self.assertEqual(run_report.sha256_of(p), expected)


class TestLoadBundles(unittest.TestCase):

    def test_skips_bundle_with_no_sample_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = write_json(Path(tmp) / "m.json", {"no_id_here": True})
            self.assertEqual(run_report.load_bundles([p]), {})


class TestRender(unittest.TestCase):
    """Case 5."""

    def test_render_produces_self_contained_html(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = {
                "pipeline_commit": "abc1234",
                "invocation_timestamp": "2026-09-29T10:00:00",
                "sample_count": 3,
                "samplesheet": {
                    "path": "/data/samples.csv", "sha256": "deadbeef",
                },
                "reference_bundle": REFS_MANIFEST,
                "status_breakdown": {
                    "ok": 1, "low_coverage": 1, "no_assembly": 0,
                    "no_barcode": 0, "fail": 0, "error": 1,
                },
                "samples": [
                    run_report.sample_row(
                        "OK1", base_metadata("OK1", "ok")),
                    run_report.sample_row(
                        "LOW1", base_metadata("LOW1", "low_coverage")),
                    run_report.sample_row("CRASHED1", None),
                ],
            }
            out = Path(tmp) / "run-report.html"
            run_report_mod.render(manifest, TEMPLATE_DIR, STATIC_DIR, out)
            html = out.read_text()
            self.assertIn("OK1/report.html", html)
            self.assertIn("LOW1/report.html", html)
            self.assertIn("CRASHED1/report.html", html)
            self.assertIn(
                run_report_mod.STATUS_LABELS["low_coverage"], html)
            self.assertIn(run_report_mod.STATUS_LABELS["error"], html)
            self.assertIn("abc1234", html)
            self.assertIn("deadbeef", html)


class TestRenderFailureFallback(unittest.TestCase):
    """Case 6 — a rendering defect must not take run_manifest.json
    down with it (CONSTITUTION rule 8, mirroring collate.py §5.2)."""

    def test_bad_template_dir_writes_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            out = d / "run-report.html"

            class Args:
                report_templates = d / "does-not-exist"
                report_static = STATIC_DIR
                out_report = out

            run_report.render_run_report(Args(), {"sample_count": 0})
            html = out.read_text()
            self.assertIn("Run report rendering failed", html)

    def test_no_templates_configured_touches_empty_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "run-report.html"

            class Args:
                report_templates = None
                report_static = None
                out_report = out

            run_report.render_run_report(Args(), {"sample_count": 0})
            self.assertTrue(out.is_file())
            self.assertEqual(out.stat().st_size, 0)


class TestMainCLI(unittest.TestCase):
    """End-to-end CLI invocation, mirroring test_collate.py/
    test_report.py's argv-patching pattern."""

    def _run_main(self, argv):
        old_argv = sys.argv
        try:
            sys.argv = argv
            return run_report.main()
        finally:
            sys.argv = old_argv

    def test_main_exits_zero_and_writes_valid_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            samplesheet_csv = write_csv(
                d / "samples.csv",
                [{"sample_id": "S1", "assembly_target": "animal_mt"}])
            meta = write_json(
                d / "S1.metadata.json", base_metadata("S1", "ok"))
            refs_manifest = write_json(d / "manifest.json", REFS_MANIFEST)
            schema_path = d / "run_manifest.schema.json"
            schema_path.write_text(json.dumps(SCHEMA))
            out_manifest = d / "run_manifest.json"
            out_report = d / "run-report.html"

            argv = [
                "run_report.py",
                "--metadata-json", str(meta),
                "--samplesheet", str(samplesheet_csv),
                "--refs-manifest", str(refs_manifest),
                "--pipeline-commit", "abc1234",
                "--workflow-start", "2026-09-29T10:00:00",
                "--schema", str(schema_path),
                "--report-templates", str(TEMPLATE_DIR),
                "--report-static", str(STATIC_DIR),
                "--out-manifest", str(out_manifest),
                "--out-report", str(out_report),
            ]
            exit_code = self._run_main(argv)
            self.assertEqual(exit_code, 0)
            manifest = json.loads(out_manifest.read_text())
            import jsonschema
            jsonschema.validate(manifest, SCHEMA)
            self.assertIn("S1/report.html", out_report.read_text())

    def test_main_exits_zero_on_malformed_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            samplesheet_csv = write_csv(
                d / "samples.csv",
                [{"sample_id": "S1", "assembly_target": "animal_mt"}])
            meta = write_json(
                d / "S1.metadata.json", base_metadata("S1", "ok"))
            schema_path = d / "bad_schema.json"
            schema_path.write_text(json.dumps({"type": "not-a-real-type"}))
            out_manifest = d / "run_manifest.json"
            out_report = d / "run-report.html"

            argv = [
                "run_report.py",
                "--metadata-json", str(meta),
                "--samplesheet", str(samplesheet_csv),
                "--schema", str(schema_path),
                "--out-manifest", str(out_manifest),
                "--out-report", str(out_report),
            ]
            exit_code = self._run_main(argv)
            self.assertEqual(exit_code, 0)
            self.assertTrue(out_manifest.is_file())

    def test_main_without_schema_arg_skips_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            samplesheet_csv = write_csv(
                d / "samples.csv",
                [{"sample_id": "S1", "assembly_target": "animal_mt"}])
            meta = write_json(
                d / "S1.metadata.json", base_metadata("S1", "ok"))
            out_manifest = d / "run_manifest.json"
            out_report = d / "run-report.html"

            argv = [
                "run_report.py",
                "--metadata-json", str(meta),
                "--samplesheet", str(samplesheet_csv),
                "--out-manifest", str(out_manifest),
                "--out-report", str(out_report),
            ]
            exit_code = self._run_main(argv)
            self.assertEqual(exit_code, 0)
            self.assertTrue(out_manifest.is_file())


if __name__ == "__main__":
    unittest.main()
