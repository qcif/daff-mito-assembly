"""Unit tests for the report renderer (task 43a) — `bin/report/` (the
Jinja machinery, `render()`, `key_findings()`, `build_warnings()`) and
the `read_qc` / recruitment plumbing added to `bin/collate.py`.

Real templates and a real Jinja render, not mocks (task 27 §6's
boundary-mocking rule only bites at an external tool call, and this
module has none) — fixtures are metadata.json-shaped dicts mirroring
the real schema, rendered through the actual `scripts/report/`
templates/static tree.

Cases (task 43a §9, numbered to match):
  1. All five sample_status values render without raising, each
     producing a Key findings block naming its own state.
  2. no_assembly and no_barcode produce different text.
  3. low_coverage Key findings does not use failure wording and is
     distinguishable from ok.
  4. Self-containment: no src="http, href="http or url(http outside
     the one deliberate informational link (the repo subtitle).
  5. read_qc parsing: well-formed, malformed and absent NanoStats.txt;
     filter_yield maths including the zero-reads edge case.
  6. Missing optional artefacts (zero-byte SVG, absent PNG, absent
     GFF, absent read_qc) render without raising.
  7. Render failure produces the fallback report and collate.py still
     exits 0.
  8. Emitted metadata.json still validates against the extended
     schema with read_qc/recruitment populated.
"""

import importlib.util
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
REPO_ROOT = BIN_DIR.parent
TEMPLATE_DIR = REPO_ROOT / "scripts" / "report" / "templates"
STATIC_DIR = REPO_ROOT / "scripts" / "report" / "static"

sys.path.insert(0, str(BIN_DIR))

import report.report as report_mod  # noqa: E402

SCHEMA = json.loads(
    (REPO_ROOT / "assets" / "sample_metadata.schema.json").read_text())

_spec = importlib.util.spec_from_file_location(
    "collate", BIN_DIR / "collate.py")
collate = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(collate)


def base_metadata(status: str, **overrides) -> dict:
    metadata = {
        "$schema": "wf5/sample-metadata/v1",
        "sample_id": "TESTSAMPLE01",
        "assembly_target": "animal_mt",
        "kingdom": "animal",
        "organelle": "mt",
        "sample_info": "",
        "sample_type": "",
        "sample_receipt_date": "",
        "storage_location": "",
        "sample_status": status,
        "sample_status_reason": "unit test reason",
        "bundle": "full" if status != "fail" else "minimal",
        "coverage": {
            "gate": {
                "estimated_cov": 42.0,
                "status": "ok",
                "coverage_basis": "target_assigned",
            },
            "estimate": {"post_subsample_cov": 42.0},
            "recruitment": {"reads_aligned": 100, "reads_recruited": 90},
        },
        "assembly": {
            "contig_count": 1, "n50": 16000, "total_bp": 16000,
            "contigs": [], "bin_metadata": None,
        },
        "homology": {
            "top_hits": [{
                "qaccver": "contig_1", "saccver": "NC_000001",
                "pident": 98.4, "length": 16000, "qcovs": 99.0,
                "evalue": 0.0, "bitscore": 5000.0,
                "stitle": "Test species mitochondrion",
            }],
        },
        "barcodes": {
            "loci": [{"locus": "COX1", "status": "pass"}],
            "n_passed": 1,
        },
        "annotation": {
            "status": "ok",
            "feature_counts": {"gene": 37, "CDS": 13, "tRNA": 22, "rRNA": 2},
            "genetic_code_agreement": True,
            "cds_crosscheck": {"agreed": [], "annotator_only": []},
        },
        "read_qc": None,
        "provenance": {
            "pipeline_commit": "abc1234",
            "reference_bundle_version": "v2026.08_1",
            "reference_bundle_generated_at": "2026-08-01T00:00:00Z",
            "tool_versions": {"flye": "2.9.6"},
            "genetic_code": {"table": 5},
        },
    }
    metadata.update(overrides)
    return metadata


class TestKeyFindingsAllStatuses(unittest.TestCase):

    def test_all_five_statuses_render_and_produce_findings(self):
        for status in (
            "ok", "low_coverage", "no_assembly", "no_barcode", "fail",
        ):
            metadata = base_metadata(status)
            findings = report_mod.key_findings(metadata)
            self.assertTrue(findings, msg=status)
            for f in findings:
                self.assertIn("text", f)
                self.assertIn("class", f)

    def test_no_assembly_and_no_barcode_produce_different_text(self):
        no_assembly = report_mod.key_findings(base_metadata("no_assembly"))
        no_barcode = report_mod.key_findings(base_metadata("no_barcode"))
        self.assertNotEqual(
            [f["text"] for f in no_assembly],
            [f["text"] for f in no_barcode],
        )
        # no_assembly must not claim an assembly outcome/annotation —
        # it is a single terminal finding, unlike no_barcode's several.
        self.assertEqual(len(no_assembly), 1)
        self.assertGreater(len(no_barcode), 1)

    def test_low_coverage_is_not_framed_as_failure(self):
        low_coverage = report_mod.key_findings(base_metadata("low_coverage"))
        ok = report_mod.key_findings(base_metadata("ok"))
        low_text = " ".join(f["text"] for f in low_coverage).lower()
        low_classes = {f["class"] for f in low_coverage}
        self.assertNotIn("danger", low_classes)
        self.assertNotIn("no organelle was assembled", low_text)
        self.assertIn("partial result", low_text)
        self.assertNotEqual(
            [f["text"] for f in low_coverage], [f["text"] for f in ok],
        )

    def test_fail_and_no_assembly_are_danger_class(self):
        for status in ("fail", "no_assembly"):
            findings = report_mod.key_findings(base_metadata(status))
            self.assertEqual(findings[0]["class"], "danger")


class TestBuildWarnings(unittest.TestCase):

    def test_ok_sample_has_no_warnings(self):
        self.assertEqual(report_mod.build_warnings(base_metadata("ok")), [])

    def test_low_coverage_warns_on_its_own_tab(self):
        warnings = report_mod.build_warnings(base_metadata("low_coverage"))
        self.assertTrue(
            any(w["tab"] == "Recruitment & coverage gate" for w in warnings)
        )

    def test_annotator_failed_warns(self):
        metadata = base_metadata("ok", annotation={
            "status": "annotator_failed",
            "reason": "mitos2 exited 1",
            "feature_counts": {"gene": 13},
        })
        warnings = report_mod.build_warnings(metadata)
        self.assertTrue(any(w["severity"] == "danger" for w in warnings))

    def test_genetic_code_disagreement_warns(self):
        metadata = base_metadata("ok")
        metadata["annotation"]["genetic_code_agreement"] = False
        warnings = report_mod.build_warnings(metadata)
        self.assertTrue(
            any("genetic-code" in w["text"] for w in warnings)
        )

    def test_total_recruited_basis_warns(self):
        metadata = base_metadata("ok")
        metadata["coverage"]["gate"]["coverage_basis"] = "total_recruited"
        warnings = report_mod.build_warnings(metadata)
        self.assertTrue(
            any(w["tab"] == "Recruitment & coverage gate" for w in warnings)
        )

    def test_cds_crosscheck_disagreement_warns(self):
        metadata = base_metadata("ok")
        metadata["annotation"]["cds_crosscheck"] = {
            "annotator_only": [{"gene": "nad1", "reason": "no miniprot hit"}],
        }
        warnings = report_mod.build_warnings(metadata)
        self.assertTrue(any("disagree" in w["text"] for w in warnings))

    def test_no_barcode_warns_on_its_own_tab(self):
        warnings = report_mod.build_warnings(base_metadata("no_barcode"))
        self.assertTrue(
            any(w["tab"] == "Extracted barcode panel" for w in warnings)
        )


class TestKeyFindingsSparseMetadata(unittest.TestCase):
    """Every helper's "no data available" branch — a sample that
    cleared the gate but carries none of the optional detail yet
    (assembly/homology/annotation/barcodes all null or empty)."""

    def _sparse(self, status):
        return {
            "sample_id": "SPARSE01",
            "sample_status": status,
            "sample_status_reason": "reason",
            "coverage": {"gate": {}, "estimate": {}},
            "assembly": {},
            "homology": {},
            "barcodes": {},
            "annotation": {},
        }

    def test_ok_with_no_optional_data(self):
        findings = report_mod.key_findings(self._sparse("ok"))
        texts = " ".join(f["text"] for f in findings)
        self.assertIn("An organelle assembly was produced.", texts)
        self.assertIn("Coverage estimate unavailable.", texts)
        self.assertIn("No annotation was produced.", texts)
        self.assertIn("No barcode loci were evaluated.", texts)
        self.assertNotIn("Top BLAST hit", texts)

    def test_low_coverage_with_no_optional_data(self):
        findings = report_mod.key_findings(self._sparse("low_coverage"))
        self.assertTrue(findings)

    def test_no_barcode_with_no_optional_data(self):
        findings = report_mod.key_findings(self._sparse("no_barcode"))
        self.assertTrue(findings)


class TestRenderSelfContainment(unittest.TestCase):

    # The repo's own subtitle link, plus a "Produced with Plotly" link
    # baked as a literal string into the vendored, inert plotly-basic
    # JS source (not a network fetch — no chart calls it at render
    # time yet). Neither is something this task's edits could add or
    # remove; both predate it.
    ALLOWED_EXTERNAL = (
        "github.com/qcif/daff-biosecurity-wf5",
        "https://plotly.com/",
    )

    def _render(self, metadata, **kwargs):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "report.html"
            report_mod.render(
                metadata, TEMPLATE_DIR, STATIC_DIR, out, **kwargs)
            return out.read_text()

    def test_self_contained_for_every_status(self):
        pattern = re.compile(r'(?:src|href)="http[^"]*"|url\(http[^)]*\)')
        for status in (
            "ok", "low_coverage", "no_assembly", "no_barcode", "fail",
        ):
            html = self._render(base_metadata(status))
            externals = [
                m for m in pattern.findall(html)
                if not any(allowed in m for allowed in self.ALLOWED_EXTERNAL)
            ]
            self.assertEqual(externals, [], msg=status)

    def test_sample_id_and_status_appear(self):
        html = self._render(base_metadata("ok"))
        self.assertIn("TESTSAMPLE01", html)

    def test_file_src_helpers_handle_absence_and_zero_byte(self):
        with tempfile.TemporaryDirectory() as tmp:
            zero = Path(tmp) / "zero.png"
            zero.write_bytes(b"")
            self.assertIsNone(report_mod._file_src(None, "image/png"))
            self.assertIsNone(report_mod._file_src(zero, "image/png"))
            self.assertIsNone(
                report_mod._file_src(Path(tmp) / "absent.png", "image/png"))

            nanoplot_dir = Path(tmp) / "np"
            nanoplot_dir.mkdir()
            self.assertIsNone(report_mod._nanoplot_report_src(None))
            self.assertIsNone(
                report_mod._nanoplot_report_src(nanoplot_dir))
            (nanoplot_dir / "NanoPlot-report.html").write_bytes(b"")
            self.assertIsNone(
                report_mod._nanoplot_report_src(nanoplot_dir))

            png = Path(tmp) / "graph.png"
            png.write_bytes(b"\x89PNG\r\n")
            self.assertTrue(
                report_mod._file_src(png, "image/png")
                .startswith("data:image/png;base64,"))

    def test_missing_optional_artefacts_do_not_raise(self):
        with tempfile.TemporaryDirectory() as tmp:
            zero_svg = Path(tmp) / "map.svg"
            zero_svg.write_bytes(b"")
            html = self._render(
                base_metadata("ok"),
                organelle_map_svg=zero_svg,
                graph_png=Path(tmp) / "absent.png",
                annotation_gff=Path(tmp) / "absent.gff",
                nanoplot_raw_dir=None,
                nanoplot_clean_dir=None,
            )
            self.assertIn("TESTSAMPLE01", html)

    def test_organelle_map_svg_wired_into_context(self):
        # 43a wires the artefact into the render context (§4.3); 43b
        # is what puts it on a tab pane. Assert on the context, not on
        # visible HTML output that doesn't exist yet.
        with tempfile.TemporaryDirectory() as tmp:
            svg = Path(tmp) / "map.svg"
            svg.write_text('<svg><circle id="gene1"/></svg>')
            context = report_mod.build_context(
                base_metadata("ok"), params={},
                nanoplot_raw_dir=None, nanoplot_clean_dir=None,
                organelle_map_svg=svg, graph_png=None, annotation_gff=None,
            )
            self.assertEqual(
                context["organelle_map_svg"],
                '<svg><circle id="gene1"/></svg>',
            )

    def test_nanoplot_report_wired_as_data_uri_in_context(self):
        with tempfile.TemporaryDirectory() as tmp:
            nanoplot_dir = Path(tmp) / "nanoplot_raw"
            nanoplot_dir.mkdir()
            (nanoplot_dir / "NanoPlot-report.html").write_text("<html></html>")
            context = report_mod.build_context(
                base_metadata("ok"), params={},
                nanoplot_raw_dir=nanoplot_dir, nanoplot_clean_dir=None,
                organelle_map_svg=None, graph_png=None, annotation_gff=None,
            )
            self.assertTrue(
                context["nanoplot_raw_src"].startswith(
                    "data:text/html;base64,"))


class TestReadQCParsing(unittest.TestCase):

    NANOSTATS = (
        "number_of_reads\t1242\n"
        "number_of_bases\t7900845.0\n"
        "median_read_length\t3374.0\n"
        "mean_read_length\t6361.4\n"
        "read_length_stdev\t7065.2\n"
        "n50\t13081.0\n"
        "mean_qual\t11.0\n"
        "median_qual\t11.5\n"
        "Reads >Q10:\t1002 (80.7%) 6.1Mb\n"
    )

    def test_well_formed_stats(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "NanoStats.txt").write_text(self.NANOSTATS)
            stats = collate.parse_nanostats(d)
            self.assertEqual(stats["number_of_reads"], 1242)
            self.assertAlmostEqual(stats["n50"], 13081.0)
            self.assertNotIn("Reads >Q10", stats)

    def test_missing_directory(self):
        self.assertIsNone(collate.parse_nanostats(None))

    def test_missing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(collate.parse_nanostats(Path(tmp)))

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "NanoStats.txt").touch()
            self.assertIsNone(collate.parse_nanostats(d))

    def test_malformed_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "NanoStats.txt").write_text(
                "number_of_reads\tnot-a-number\ngarbage line with no tab\n"
            )
            self.assertIsNone(collate.parse_nanostats(d))

    def test_filter_yield_maths(self):
        raw = {"number_of_reads": 1000, "number_of_bases": 5_000_000}
        clean = {"number_of_reads": 800, "number_of_bases": 4_000_000}
        yield_ = collate.filter_yield(raw, clean)
        self.assertEqual(yield_["reads_retained"], 800)
        self.assertEqual(yield_["reads_retained_pct"], 80.0)
        self.assertEqual(yield_["bases_retained_pct"], 80.0)

    def test_filter_yield_zero_reads_edge_case(self):
        raw = {"number_of_reads": 0, "number_of_bases": 0}
        clean = {"number_of_reads": 0, "number_of_bases": 0}
        yield_ = collate.filter_yield(raw, clean)
        self.assertIsNone(yield_["reads_retained_pct"])
        self.assertIsNone(yield_["bases_retained_pct"])

    def test_filter_yield_missing_side_returns_none(self):
        self.assertIsNone(collate.filter_yield(None, {"number_of_reads": 1}))
        self.assertIsNone(collate.filter_yield({"number_of_reads": 1}, None))

    def test_read_qc_section_absent_both(self):
        self.assertIsNone(collate.read_qc_section(None, None))

    def test_read_qc_section_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw_dir = Path(tmp) / "raw"
            clean_dir = Path(tmp) / "clean"
            raw_dir.mkdir()
            clean_dir.mkdir()
            (raw_dir / "NanoStats.txt").write_text(self.NANOSTATS)
            (clean_dir / "NanoStats.txt").write_text(self.NANOSTATS)
            section = collate.read_qc_section(raw_dir, clean_dir)
            self.assertIsNotNone(section["raw"])
            self.assertIsNotNone(section["clean"])
            self.assertIsNotNone(section["filter_yield"])


class TestRecruitmentSection(unittest.TestCase):

    def test_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "recruit_stats.json"
            path.write_text(json.dumps({"reads_aligned": 10}))
            self.assertEqual(
                collate.recruitment_section(path), {"reads_aligned": 10})

    def test_absent(self):
        self.assertIsNone(collate.recruitment_section(None))


class TestRenderFailureFallback(unittest.TestCase):
    """§5.2 / §9 case 7 — a rendering defect must not take the sample's
    bundle down; COLLATE must still exit 0 with a fallback report."""

    def test_render_bundle_report_writes_fallback_on_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "report.html"

            class Args:
                report_templates = Path(tmp) / "does-not-exist"
                report_static = STATIC_DIR
                params_json = None
                nanoplot_raw = None
                nanoplot_clean = None
                organelle_map_svg = None
                graph_png = None
                annotation_gff = None
                out_report = out

            collate.render_bundle_report(Args(), base_metadata("ok"))
            html = out.read_text()
            self.assertIn("Report rendering failed", html)
            self.assertIn("TESTSAMPLE01", html)

    def test_no_templates_configured_touches_empty_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "report.html"

            class Args:
                report_templates = None
                report_static = None
                out_report = out

            collate.render_bundle_report(Args(), base_metadata("ok"))
            self.assertTrue(out.is_file())
            self.assertEqual(out.stat().st_size, 0)

    def test_collate_main_exits_zero_on_render_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            status_json = d / "status.json"
            status_json.write_text(json.dumps({"status": "ok"}))
            coverage_json = d / "coverage.json"
            coverage_json.write_text(json.dumps({}))
            meta_json = d / "meta.json"
            meta_json.write_text(json.dumps({
                "sample_id": "S1", "assembly_target": "animal_mt",
            }))
            bad_templates = d / "not-a-real-template-dir"
            out_report = d / "report.html"
            out_metadata = d / "metadata.json"

            argv = [
                "collate.py",
                "--meta-json", str(meta_json),
                "--status-json", str(status_json),
                "--coverage-json", str(coverage_json),
                "--report-templates", str(bad_templates),
                "--report-static", str(STATIC_DIR),
                "--out-metadata", str(out_metadata),
                "--out-report", str(out_report),
            ]
            old_argv = sys.argv
            old_cwd = Path.cwd()
            try:
                sys.argv = argv
                import os
                os.chdir(tmp)
                rc = collate.main()
            finally:
                sys.argv = old_argv
                import os
                os.chdir(old_cwd)
            self.assertEqual(rc, 0)
            self.assertIn("Report rendering failed", out_report.read_text())


class TestCssHashFilter(unittest.TestCase):
    """`css_hash` is registered as a Jinja filter (task 43a §2 — kept
    from the boilerplate for 43b's use) but no template calls it yet;
    direct unit coverage per rule 14."""

    def test_deterministic_and_length(self):
        from report.filters.css_hash import css_hash
        self.assertEqual(css_hash("same-input"), css_hash("same-input"))
        self.assertEqual(len(css_hash("x")), 10)
        self.assertEqual(len(css_hash("x", length=6)), 6)


class TestSchemaValidation(unittest.TestCase):

    def test_metadata_with_read_qc_and_recruitment_validates(self):
        import jsonschema
        metadata = base_metadata("ok", read_qc={
            "raw": {"number_of_reads": 100},
            "clean": {"number_of_reads": 90},
            "filter_yield": {"reads_retained": 90},
        })
        jsonschema.validate(metadata, SCHEMA)

    def test_metadata_with_null_read_qc_validates(self):
        import jsonschema
        jsonschema.validate(base_metadata("ok"), SCHEMA)


if __name__ == "__main__":
    unittest.main()
