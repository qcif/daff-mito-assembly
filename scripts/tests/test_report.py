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

import report.confidence as confidence  # noqa: E402
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
            "loci": [{
                "gene": "COX1", "status": "pass", "reason": "",
                "seqid": "contig_1", "start": "1", "end": "600",
                "strand": "+", "identity": "0.95", "genetic_code": "5",
                "length_nt": "600",
            }],
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
            any(w["tab"] == "Validation" for w in warnings)
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
            any(w["tab"] == "Validation" for w in warnings)
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
            any(w["tab"] == "Barcodes" for w in warnings)
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
                barcodes_fasta=None, workflow_start=None,
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
                barcodes_fasta=None, workflow_start=None,
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
                barcodes_fasta = None
                workflow_start = None
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


# ---------------------------------------------------------------------------
# Task 43b — four-tab restructure. Cases numbered against task
# 43b_report_stage_tabs.md §6.
# ---------------------------------------------------------------------------


INT_ANIMAL_METADATA = json.loads(
    (REPO_ROOT / "tests" / "integration" / "output" / "INT-ANIMAL-01"
     / "metadata.json").read_text()
)
INT_PLANT_PT_METADATA = json.loads(
    (REPO_ROOT / "tests" / "integration" / "output" / "INT-PLANT-01-pt"
     / "metadata.json").read_text()
)


class TestOrganelleTitle(unittest.TestCase):

    def test_mt_title(self):
        context = report_mod.build_context(
            base_metadata("ok", organelle="mt"), params={},
            nanoplot_raw_dir=None, nanoplot_clean_dir=None,
            organelle_map_svg=None, graph_png=None, annotation_gff=None,
            barcodes_fasta=None, workflow_start=None,
        )
        self.assertEqual(context["title"], "Mitochondrial genome assembly")

    def test_pt_title(self):
        context = report_mod.build_context(
            base_metadata("ok", organelle="pt"), params={},
            nanoplot_raw_dir=None, nanoplot_clean_dir=None,
            organelle_map_svg=None, graph_png=None, annotation_gff=None,
            barcodes_fasta=None, workflow_start=None,
        )
        self.assertEqual(context["title"], "Chloroplast genome assembly")

    def test_unknown_organelle_falls_back_to_default(self):
        context = report_mod.build_context(
            base_metadata("ok", organelle="unknown"), params={},
            nanoplot_raw_dir=None, nanoplot_clean_dir=None,
            organelle_map_svg=None, graph_png=None, annotation_gff=None,
            barcodes_fasta=None, workflow_start=None,
        )
        self.assertEqual(context["title"], report_mod.DEFAULT_TITLE)
        self.assertNotIn("Mitochondrial", context["title"])
        self.assertNotIn("Chloroplast", context["title"])


class TestWallTime(unittest.TestCase):

    def test_absent_start_renders_none(self):
        ctx = report_mod.wall_time_context(None)
        self.assertIsNone(ctx["start_time"])
        self.assertIsNone(ctx["end_time"])
        self.assertIsNone(ctx["duration"])

    def test_malformed_start_falls_back(self):
        ctx = report_mod.wall_time_context("not-a-timestamp")
        self.assertEqual(ctx["start_time"], "not-a-timestamp")
        self.assertIsNone(ctx["duration"])

    def test_valid_start_computes_duration(self):
        from datetime import datetime, timedelta
        start = (datetime.now() - timedelta(hours=1, minutes=2, seconds=3))
        ctx = report_mod.wall_time_context(start.isoformat())
        self.assertIsNotNone(ctx["start_time"])
        self.assertIsNotNone(ctx["end_time"])
        self.assertRegex(ctx["duration"], r"^\d{2}:\d{2}:\d{2}$")

    def test_facility_analyst_unset_render_dash(self):
        context = report_mod.build_context(
            base_metadata("ok"), params={},
            nanoplot_raw_dir=None, nanoplot_clean_dir=None,
            organelle_map_svg=None, graph_png=None, annotation_gff=None,
            barcodes_fasta=None, workflow_start=None,
        )
        self.assertEqual(context["facility"], "-")
        self.assertEqual(context["analyst_name"], "-")

    def test_facility_analyst_set_from_params(self):
        context = report_mod.build_context(
            base_metadata("ok"),
            params={"facility": "QCIF Lab", "analyst_name": "A. Analyst"},
            nanoplot_raw_dir=None, nanoplot_clean_dir=None,
            organelle_map_svg=None, graph_png=None, annotation_gff=None,
            barcodes_fasta=None, workflow_start=None,
        )
        self.assertEqual(context["facility"], "QCIF Lab")
        self.assertEqual(context["analyst_name"], "A. Analyst")

    def test_facility_analyst_dash_rendered_in_html(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "report.html"
            report_mod.render(
                base_metadata("ok"), TEMPLATE_DIR, STATIC_DIR, out,
                params={},
            )
            html = out.read_text()
            self.assertIn(">-</td>", html)


class TestConfidenceClassifier(unittest.TestCase):
    """Direct unit coverage of the quadrant classifier, independent of
    rendering (task 43b §6 — "the quadrant classifier is unit-tested
    directly, at and either side of each threshold")."""

    def test_null_triplet_is_unscored(self):
        self.assertEqual(
            confidence.classify(None, 90.0), confidence.UNSCORED)
        self.assertEqual(
            confidence.classify(80.0, None), confidence.UNSCORED)
        self.assertEqual(
            confidence.classify(None, None), confidence.UNSCORED)

    def test_complete_well_referenced(self):
        self.assertEqual(
            confidence.classify(80.0, 90.0),
            confidence.COMPLETE_WELL_REFERENCED)

    def test_truncated_well_referenced(self):
        self.assertEqual(
            confidence.classify(68.75, 13.0),
            confidence.TRUNCATED_WELL_REFERENCED)

    def test_complete_under_referenced_real_examples(self):
        # INT-ANIMAL-01's real spread (task 43b §5.4) — COX3 and ATP6
        # are intact genes in an under-referenced clade, not failures.
        self.assertEqual(
            confidence.classify(39.453, 96.0),
            confidence.COMPLETE_UNDER_REFERENCED)
        self.assertEqual(
            confidence.classify(39.286, 100.0),
            confidence.COMPLETE_UNDER_REFERENCED)

    def test_fragmentary(self):
        self.assertEqual(
            confidence.classify(38.863, 68.0), confidence.FRAGMENTARY)

    def test_threshold_boundaries_pident(self):
        below = confidence.PIDENT_THRESHOLD - 0.1
        at = confidence.PIDENT_THRESHOLD
        above = confidence.PIDENT_THRESHOLD + 0.1
        qcov = confidence.QCOVHSP_THRESHOLD + 10
        self.assertEqual(
            confidence.classify(below, qcov),
            confidence.COMPLETE_UNDER_REFERENCED)
        self.assertEqual(
            confidence.classify(at, qcov),
            confidence.COMPLETE_WELL_REFERENCED)
        self.assertEqual(
            confidence.classify(above, qcov),
            confidence.COMPLETE_WELL_REFERENCED)

    def test_threshold_boundaries_qcovhsp(self):
        pident = confidence.PIDENT_THRESHOLD + 10
        below = confidence.QCOVHSP_THRESHOLD - 0.1
        at = confidence.QCOVHSP_THRESHOLD
        above = confidence.QCOVHSP_THRESHOLD + 0.1
        self.assertEqual(
            confidence.classify(pident, below),
            confidence.TRUNCATED_WELL_REFERENCED)
        self.assertEqual(
            confidence.classify(pident, at),
            confidence.COMPLETE_WELL_REFERENCED)
        self.assertEqual(
            confidence.classify(pident, above),
            confidence.COMPLETE_WELL_REFERENCED)

    def test_label_and_description_cover_every_quadrant(self):
        for q in (
            confidence.COMPLETE_WELL_REFERENCED,
            confidence.TRUNCATED_WELL_REFERENCED,
            confidence.COMPLETE_UNDER_REFERENCED,
            confidence.FRAGMENTARY,
            confidence.UNSCORED,
        ):
            self.assertTrue(confidence.label(q))
            self.assertTrue(confidence.description(q))

    def test_label_unknown_quadrant_falls_back_to_value(self):
        self.assertEqual(confidence.label("bogus"), "bogus")
        self.assertEqual(confidence.description("bogus"), "")


class TestNoCompositeScore(unittest.TestCase):
    """§5.4 rule 1 — the regression that matters most: no composite
    score exists anywhere. Mechanically checkable, not merely a style
    guideline."""

    def test_score_row_never_combines_axes(self):
        row = report_mod._score_row({
            "gene": "COX3", "pident": 39.453, "qcovhsp": 96.0,
            "bitscore": 196.0, "seqid": "contig_10", "start": 419,
        })
        # bitscore/pident/qcovhsp are legitimate single-axis fields
        # from upstream tools, not composites of the two confidence
        # axes — excluded from the forbidden-substring scan below.
        allowed = {"bitscore", "pident", "qcovhsp"}
        forbidden = (
            "confidence", "score", "combined", "mean", "weighted",
            "average",
        )
        for key in row:
            if key in allowed:
                continue
            lowered = key.lower()
            for word in forbidden:
                self.assertNotIn(
                    word, lowered,
                    msg=f"result row exposes a composite-looking key: {key}")
        # pident and qcovhsp are exposed as two separately-labelled
        # values, not reduced to one.
        self.assertIn("pident", row)
        self.assertIn("qcovhsp", row)
        self.assertNotEqual(row["pident"], row["qcovhsp"])

    def test_rendered_row_shows_both_axes_separately_labelled(self):
        html = _render(base_metadata_with_scores())
        self.assertIn("Identity (pident)", html)
        self.assertIn("Coverage (qcovhsp)", html)

    def test_low_identity_high_coverage_carries_affirmative_label(self):
        # COX3 from INT-ANIMAL-01, 39.5% identity / 96% coverage — a
        # real intact gene that must not read as a single-scale
        # failure (task 43b §6).
        html = _render(INT_ANIMAL_METADATA)
        self.assertIn("Complete gene, under-referenced clade", html)

    def test_null_triplet_renders_unscored_not_zero(self):
        metadata = base_metadata_with_scores()
        metadata["annotation"]["cds_scores"] = [{
            "gene": "ND9", "seqid": "contig_1", "start": 1, "end": 100,
            "source": "mitos", "pident": None, "qcovhsp": None,
            "bitscore": None,
        }]
        html = _render(metadata)
        self.assertIn("unscored", html)
        # No row in this render should show a lone "0" as the
        # confidence value for the unscored feature.
        self.assertNotIn(">0%</td>", html)

    def test_bitscore_does_not_drive_cross_gene_sort_order(self):
        # A high-bitscore long gene (COX1-like) must not outrank a
        # better-covered short gene purely on bitscore — sort is by
        # genomic position only (task 43b §5.4).
        metadata = base_metadata_with_scores()
        metadata["annotation"]["cds_scores"] = [
            {
                "gene": "HIGH_BITSCORE_LONG", "seqid": "contig_1",
                "start": 500, "end": 2500, "source": "miniprot",
                "pident": 70.0, "qcovhsp": 99.0, "bitscore": 900.0,
            },
            {
                "gene": "BETTER_COVERED_SHORT", "seqid": "contig_1",
                "start": 10, "end": 110, "source": "mitos",
                "pident": 60.0, "qcovhsp": 85.0, "bitscore": 60.0,
            },
        ]
        view = report_mod.assembly_view(metadata)
        genes_in_order = [s["gene"] for s in view["cds_scores"]]
        # Sorted by (seqid, start): the short gene at start=10 comes
        # first, even though its bitscore (60.0) is far below the
        # other gene's (900.0) — bitscore never enters the sort key.
        self.assertEqual(
            genes_in_order, ["BETTER_COVERED_SHORT", "HIGH_BITSCORE_LONG"])


def base_metadata_with_scores(**overrides):
    metadata = base_metadata("ok", **overrides)
    metadata["annotation"] = {
        **metadata["annotation"],
        "cds_scores": [{
            "gene": "COX3", "seqid": "contig_1", "start": 100,
            "end": 900, "source": "miniprot",
            "pident": 39.453, "qcovhsp": 96.0, "bitscore": 196.0,
        }],
    }
    return metadata


def _render(metadata, **kwargs):
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "report.html"
        report_mod.render(metadata, TEMPLATE_DIR, STATIC_DIR, out, **kwargs)
        return out.read_text()


class TestAssemblyView(unittest.TestCase):

    def test_contig_bucket_target(self):
        self.assertEqual(
            report_mod._contig_bucket("target_candidate"), "target")

    def test_contig_bucket_secondary(self):
        self.assertEqual(
            report_mod._contig_bucket("secondary_target"), "secondary")

    def test_contig_bucket_off_target(self):
        self.assertEqual(
            report_mod._contig_bucket("off_target"), "off-target")
        self.assertEqual(
            report_mod._contig_bucket("sibling_organelle"), "off-target")
        self.assertEqual(report_mod._contig_bucket(None), "off-target")

    def test_coverage_chart_data_skips_none_coverage(self):
        contigs = [
            {"contig": "c1", "coverage": 10.0, "bucket": "target"},
            {"contig": "c2", "coverage": None, "bucket": "off-target"},
        ]
        chart = report_mod._coverage_chart_data(contigs)
        self.assertEqual(chart["x"], ["c1"])
        self.assertEqual(chart["y"], [10.0])
        self.assertEqual(chart["colors"], ["#2ca02c"])

    def test_real_fixture_bucket_join(self):
        view = report_mod.assembly_view(INT_ANIMAL_METADATA)
        buckets = {c["contig"]: c["bucket"] for c in view["contigs"]}
        self.assertEqual(buckets["contig_10"], "target")
        self.assertEqual(buckets["contig_5"], "off-target")

    def test_annotator_failed_flag(self):
        metadata = base_metadata("ok")
        metadata["annotation"] = {
            "status": "annotator_failed",
            "reason": "mitos2 exited 1",
            "annotator_exit_code": 1,
        }
        view = report_mod.assembly_view(metadata)
        self.assertTrue(view["annotator_failed"])
        self.assertFalse(view["cds_only"])

    def test_cds_only_flag_plant(self):
        metadata = base_metadata("ok", kingdom="plant")
        metadata["annotation"] = {"status": "ok_cds_only"}
        view = report_mod.assembly_view(metadata)
        self.assertTrue(view["cds_only"])
        self.assertFalse(view["annotator_failed"])

    def test_cds_only_flag_false_for_animal(self):
        metadata = base_metadata("ok", kingdom="animal")
        metadata["annotation"] = {"status": "ok_cds_only"}
        view = report_mod.assembly_view(metadata)
        self.assertFalse(view["cds_only"])

    def test_plastid_applied_renders_target_source(self):
        html = _render(INT_PLANT_PT_METADATA)
        self.assertIn("c4_plastid_path1", html)

    def test_plastid_withheld_renders_warning(self):
        metadata = json.loads(json.dumps(INT_PLANT_PT_METADATA))
        metadata["assembly"]["bin_metadata"][
            "plastid_canonicalisation"] = {
                "branch": "canonical",
                "edge_count": 3,
                "lsc_edge": "edge_3",
                "ir_edge": "edge_2",
                "ssc_edge": "edge_1",
                "substitution_applied": False,
                "substitution_withheld_reason": "no_target_selected",
                "non_canonical_reason": None,
            }
        html = _render(metadata)
        self.assertIn("canonical structure, no", html.lower())
        warnings = report_mod.build_warnings(metadata)
        self.assertTrue(
            any("no taxonomic support" in w["text"] for w in warnings))

    def test_plant_cds_only_warning_mirrors_to_overview(self):
        metadata = base_metadata("ok", kingdom="plant")
        metadata["annotation"] = {"status": "ok_cds_only"}
        warnings = report_mod.build_warnings(metadata)
        self.assertTrue(any("CDS-only" in w["text"] for w in warnings))
        html = _render(metadata)
        self.assertIn("CDS-only", html)


class TestValidationView(unittest.TestCase):

    def test_below_species_threshold_flagged(self):
        metadata = base_metadata("ok")
        metadata["homology"]["top_hits"][0]["pident"] = 85.0
        view = report_mod.validation_view(metadata)
        self.assertTrue(view["top_hits"][0]["below_species_threshold"])

    def test_above_species_threshold_not_flagged(self):
        view = report_mod.validation_view(base_metadata("ok"))
        self.assertFalse(view["top_hits"][0]["below_species_threshold"])

    def test_below_species_threshold_renders_flag_not_assignment(self):
        metadata = base_metadata("ok")
        metadata["homology"]["top_hits"][0]["pident"] = 85.0
        html = _render(metadata)
        self.assertIn("below species threshold", html)
        self.assertNotIn("taxonomic assignment:", html.lower())

    def test_sibling_split_present_for_plant_target(self):
        view = report_mod.validation_view(INT_PLANT_PT_METADATA)
        self.assertEqual(view["gate"]["target_assigned_bases"], 38120021)
        self.assertEqual(view["gate"]["sibling_assigned_bases"], 3415019)

    def test_total_recruited_basis_caveat_renders(self):
        metadata = base_metadata("ok")
        metadata["coverage"]["gate"]["coverage_basis"] = "total_recruited"
        html = _render(metadata)
        self.assertIn("may over-state", html)

    def test_subsampled_renders_fraction_seed_and_pre_post(self):
        metadata = base_metadata("ok")
        metadata["coverage"]["estimate"] = {
            "subsampled": True, "fraction": 0.5, "seed": 42,
            "pre_subsample_cov": 300.0, "post_subsample_cov": 150.0,
        }
        html = _render(metadata)
        self.assertIn("0.5", html)
        self.assertIn("42", html)
        self.assertIn("300.0", html)
        self.assertIn("150.0", html)

    def test_low_coverage_renders_warn_floor_warning_in_both_places(self):
        metadata = base_metadata("low_coverage")
        html = _render(metadata)
        # Both the Validation tab body and the Overview warnings
        # mirror must carry the warn-floor text (task 43b §6/§9) — at
        # least twice; Key findings' own softer framing may add a
        # third occurrence.
        self.assertGreaterEqual(
            html.count("below the coverage warn floor"), 2,
            msg="expected the warn-floor warning in both the "
                "Validation tab and the Overview mirror")
        overview_start = html.index('id="tab-0"')
        overview_end = html.index('id="tab-1"')
        validation_start = html.index('id="tab-2"')
        validation_end = html.index('id="tab-3"')
        self.assertIn(
            "below the coverage warn floor",
            html[overview_start:overview_end])
        self.assertIn(
            "below the coverage warn floor",
            html[validation_start:validation_end])

    def test_cds_crosscheck_and_genetic_code_disagreement_surface(self):
        metadata = base_metadata("ok")
        metadata["annotation"]["cds_crosscheck"] = {
            "agreed": {"raw": ["cox1_0"], "canonical": ["COX1"]},
            "miniprot_only": {"raw": [], "canonical": []},
            "coordinate_conflicts": [],
            "annotator_only": [
                {"gene": "atp8_1", "reason": "overlap",
                 "canonical": "ATP8"},
            ],
        }
        metadata["annotation"]["genetic_code_agreement"] = False
        metadata["annotation"]["genetic_code_annotate"] = 5
        metadata["annotation"]["genetic_code_cds"] = 2
        html = _render(metadata)
        self.assertIn("ATP8", html)
        self.assertIn("different genetic-code tables", html)

    def test_soft_fail_validation_is_terminal_content(self):
        metadata = base_metadata("fail")
        html = _render(metadata)
        self.assertIn("terminal content", html)


class TestBarcodesView(unittest.TestCase):

    def test_barcode_id_matches_extract_barcodes_convention(self):
        locus = {
            "gene": "COX1", "seqid": "contig_10", "start": "2863",
            "end": "4387",
        }
        self.assertEqual(
            report_mod._barcode_id(locus), "COX1_contig_10_2863_4387")

    def test_parse_fasta_absent_and_empty(self):
        self.assertEqual(report_mod._parse_fasta(None), {})
        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "empty.fasta"
            empty.touch()
            self.assertEqual(report_mod._parse_fasta(empty), {})

    def test_parse_fasta_no_header_lines(self):
        # Malformed FASTA with sequence lines but no leading '>' —
        # current_id never gets set, exercising the loop's "no record
        # seen yet" branch.
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "malformed.fasta"
            fasta.write_text("ACGTACGT\n")
            self.assertEqual(report_mod._parse_fasta(fasta), {})

    def test_parse_fasta_multi_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "barcodes.fasta"
            fasta.write_text(">a\nACGT\nACGT\n>b\nTTTT\n")
            seqs = report_mod._parse_fasta(fasta)
            self.assertEqual(seqs, {"a": "ACGTACGT", "b": "TTTT"})

    def test_passed_locus_carries_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "barcodes.fasta"
            fasta.write_text(">COX2_contig_10_2125_2792\nACGTACGT\n")
            metadata = base_metadata("ok")
            metadata["barcodes"]["loci"] = [{
                "gene": "COX2", "status": "pass", "reason": "",
                "seqid": "contig_10", "start": "2125", "end": "2792",
                "strand": "-", "identity": "0.6", "genetic_code": "5",
                "length_nt": "668",
            }]
            view = report_mod.barcodes_view(metadata, fasta)
            self.assertEqual(view["passed"][0]["sequence"], "ACGTACGT")

    def test_every_dropout_reason_renders(self):
        metadata = base_metadata("ok")
        metadata["barcodes"]["loci"] = [
            {"gene": "A", "status": "fail", "reason": "not_found"},
            {"gene": "B", "status": "fail", "reason": "invalid_length"},
            {"gene": "C", "status": "fail",
             "reason": "identity_below_floor"},
            {"gene": "D", "status": "fail",
             "reason": "internal_stop_codon"},
        ]
        html = _render(metadata)
        for reason in (
            "not_found", "invalid_length", "identity_below_floor",
            "internal_stop_codon",
        ):
            self.assertIn(reason, html, msg=reason)

    def test_barcodes_fasta_wired_as_download(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "barcodes.fasta"
            fasta.write_text(">COX2_contig_10_2125_2792\nACGT\n")
            html = _render(
                base_metadata("ok"), barcodes_fasta=fasta)
            self.assertIn('download="barcodes.fasta"', html)
            self.assertIn("data:text/plain;base64,", html)

    def test_annotation_gff_wired_as_download(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "annotation.gff"
            gff.write_text("##gff-version 3\n")
            html = _render(base_metadata("ok"), annotation_gff=gff)
            self.assertIn('download="organelle_annotation.gff"', html)

    def test_terminal_barcodes_tab_states_no_assembly(self):
        html = _render(base_metadata("no_assembly"))
        self.assertIn(
            "no barcode extraction was possible", html)


class TestFourTabStrip(unittest.TestCase):

    def test_exactly_four_tabs_overview_active(self):
        html = _render(base_metadata("ok"))
        self.assertEqual(html.count('role="tab"'), 4)
        for name in ("Overview", "Assembly", "Validation", "Barcodes"):
            self.assertIn(name, html)
        self.assertIn('id="tab-0-btn"', html)
        # Overview's button is the one carrying the "active" class.
        active_idx = html.index('class="nav-link active"')
        tab1_idx = html.index('id="tab-1-btn"')
        self.assertLess(active_idx, tab1_idx)

    def test_key_findings_is_inside_overview_tab_not_a_header_row(self):
        html = _render(base_metadata("ok"))
        overview_start = html.index('id="tab-0"')
        overview_end = html.index('id="tab-1"')
        key_findings_idx = html.index('id="key-findings"')
        self.assertGreater(key_findings_idx, overview_start)
        self.assertLess(key_findings_idx, overview_end)
        # And it is not part of an always-visible header above the
        # strip: the tab strip markup precedes it.
        tabs_idx = html.index('id="mainTabs"')
        self.assertLess(tabs_idx, key_findings_idx)


if __name__ == "__main__":
    unittest.main()
