"""Unit tests for bin/collate.py (C6) — spec §2 stage 15, §6a.5,
task 42.

Real temporary files, not mocks — C6 shells out to nothing (task 27
§6's boundary-mocking rule only bites at an external tool call, and
this module has none). Fixtures are small on-disk JSON/TSV/FASTA
trees mirroring the real upstream artifacts' shapes.

Cases (task 42 §8, numbered to match):
  1. Each of the five sample_status values dispatches to the right
     bundle (classify_sample + bundle_kind).
  2. An unrecognised gate status fails closed.
  3. low_coverage yields a full bundle with the warning present.
  4. no_barcode yields a full bundle.
  5. Withheld-substitution plant_pt case -> no_assembly, isoforms
     shipped as diagnostics.
  6. Minimal bundle emits no empty placeholder FASTA/GFF.
  7. Emitted metadata.json validates against the schema in every
     branch.
  8. Provenance round-trip: genetic_code.json survives into
     metadata.json.
  9. Missing/malformed/zero-byte inputs still exit 0 and still
     produce a schema-valid metadata.json.

Plus direct unit coverage of every helper for 100% branch coverage
per CONSTITUTION rule 14.
"""

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
REPO_ROOT = BIN_DIR.parent

# collate.py imports annotate_summary at module scope (sibling-module
# reuse, mirroring bin_target.py's use of intervals/plastid_canonicalise
# — rule 19), relying on Nextflow's bin/ staging at runtime.
sys.path.insert(0, str(BIN_DIR))

_spec = importlib.util.spec_from_file_location(
    "collate", BIN_DIR / "collate.py")
collate = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(collate)

SCHEMA = json.loads(
    (REPO_ROOT / "assets" / "sample_metadata.schema.json").read_text())

GENE_SETS = {
    "$schema": "wf5/gene-sets/v1",
    "version": "test",
    "sets": {
        "animal_mt": {
            "name": "metazoan_37",
            "protein_coding": ["ATP6", "ATP8", "COX1", "CYTB"],
            "protein_coding_aliases": {"CYTB": ["COB", "CYB"]},
            "rrna": ["rrnL", "rrnS"],
            "trna_count": 22,
        },
    },
}


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def write_json(path: Path, obj) -> Path:
    path.write_text(json.dumps(obj))
    return path


class TestClassifySample(unittest.TestCase):
    """§3.1 dispatch matrix — case 1/2/3/4."""

    def test_fail_gate(self):
        status, reason = collate.classify_sample("fail", [], True, 0)
        self.assertEqual(status, collate.STATUS_FAIL)
        self.assertIn("hard floor", reason)

    def test_unrecognised_gate_fails_closed(self):
        status, reason = collate.classify_sample("weird", [], True, 0)
        self.assertEqual(status, collate.STATUS_FAIL)
        self.assertIn("unrecognised", reason)

    def test_no_assembly_empty_contigs(self):
        status, _ = collate.classify_sample("ok", [], False, 0)
        self.assertEqual(status, collate.STATUS_NO_ASSEMBLY)

    def test_no_assembly_empty_fasta(self):
        status, _ = collate.classify_sample("ok", ["ctg1"], True, 0)
        self.assertEqual(status, collate.STATUS_NO_ASSEMBLY)

    def test_no_barcode(self):
        status, _ = collate.classify_sample("ok", ["ctg1"], False, 0)
        self.assertEqual(status, collate.STATUS_NO_BARCODE)

    def test_low_coverage_full_bundle(self):
        status, reason = collate.classify_sample(
            "low_coverage", ["ctg1"], False, 1)
        self.assertEqual(status, collate.STATUS_LOW_COVERAGE)
        self.assertIn("warn floor", reason)
        self.assertEqual(collate.bundle_kind(status), collate.BUNDLE_FULL)

    def test_ok(self):
        status, _ = collate.classify_sample("ok", ["ctg1"], False, 1)
        self.assertEqual(status, collate.STATUS_OK)
        self.assertEqual(collate.bundle_kind(status), collate.BUNDLE_FULL)

    def test_no_barcode_is_full_bundle(self):
        status, _ = collate.classify_sample("ok", ["ctg1"], False, 0)
        self.assertEqual(collate.bundle_kind(status), collate.BUNDLE_FULL)

    def test_fail_and_no_assembly_are_minimal(self):
        self.assertEqual(
            collate.bundle_kind(collate.STATUS_FAIL), collate.BUNDLE_MINIMAL)
        self.assertEqual(
            collate.bundle_kind(collate.STATUS_NO_ASSEMBLY),
            collate.BUNDLE_MINIMAL)


class TestKingdomOrganelle(unittest.TestCase):
    def test_plant_pt(self):
        self.assertEqual(
            collate.kingdom_organelle("plant_pt"), ("plant", "pt"))

    def test_animal_mt(self):
        self.assertEqual(
            collate.kingdom_organelle("animal_mt"), ("animal", "mt"))


class TestReadJson(unittest.TestCase):
    def test_none_path(self):
        self.assertIsNone(collate.read_json(None))

    def test_missing_file(self):
        self.assertIsNone(collate.read_json("/no/such/file.json"))

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(Path(td) / "x.json", "")
            self.assertIsNone(collate.read_json(p))

    def test_malformed_json(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(Path(td) / "x.json", "{not json")
            self.assertIsNone(collate.read_json(p))

    def test_valid_json(self):
        with tempfile.TemporaryDirectory() as td:
            p = write_json(Path(td) / "x.json", {"a": 1})
            self.assertEqual(collate.read_json(p), {"a": 1})


class TestReadTsv(unittest.TestCase):
    def test_none_path(self):
        self.assertEqual(collate.read_tsv(None), [])

    def test_missing_file(self):
        self.assertEqual(collate.read_tsv("/no/such/file.tsv"), [])

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(Path(td) / "x.tsv", "")
            self.assertEqual(collate.read_tsv(p), [])

    def test_rows(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(Path(td) / "x.tsv", "a\tb\n1\t2\n3\t4\n")
            self.assertEqual(
                collate.read_tsv(p),
                [{"a": "1", "b": "2"}, {"a": "3", "b": "4"}])


class TestIsEmptyFasta(unittest.TestCase):
    def test_none(self):
        self.assertTrue(collate.is_empty_fasta(None))

    def test_missing(self):
        self.assertTrue(collate.is_empty_fasta("/no/such/file.fasta"))

    def test_zero_byte(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(Path(td) / "t.fasta", "")
            self.assertTrue(collate.is_empty_fasta(p))

    def test_non_empty(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(Path(td) / "t.fasta", ">a\nACGT\n")
            self.assertFalse(collate.is_empty_fasta(p))


class TestParseAssemblyInfo(unittest.TestCase):
    def test_missing(self):
        self.assertEqual(collate.parse_assembly_info(None), [])
        self.assertEqual(collate.parse_assembly_info("/nope"), [])

    def test_rows(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(
                Path(td) / "assembly_info.txt",
                "#seq_name\tlength\tcov.\tcirc.\n"
                "ctg1\t100\t20.0\tY\n"
                "ctg2\t50\t5.0\tN\n"
                "short\t1\n",
            )
            rows = collate.parse_assembly_info(p)
            self.assertEqual(len(rows), 2)
            self.assertTrue(rows[0]["circular"])
            self.assertFalse(rows[1]["circular"])


class TestN50(unittest.TestCase):
    def test_empty(self):
        self.assertIsNone(collate.n50([]))

    def test_typical(self):
        self.assertEqual(collate.n50([100, 90, 10]), 100)

    def test_multiple_iterations_before_reaching_half(self):
        self.assertEqual(collate.n50([10, 10, 10, 10]), 10)


class TestAssemblySection(unittest.TestCase):
    def test_nothing_present(self):
        self.assertIsNone(collate.assembly_section([], None))

    def test_present(self):
        rows = [{"contig": "c1", "length": 10, "coverage": 1.0,
                 "circular": True}]
        section = collate.assembly_section(
            rows, {"contigs_selected": ["c1"]})
        self.assertEqual(section["contig_count"], 1)
        self.assertEqual(section["total_bp"], 10)


class TestTopBlastHits(unittest.TestCase):
    def test_missing(self):
        self.assertEqual(collate.top_blast_hits(None), [])
        self.assertEqual(collate.top_blast_hits("/nope"), [])

    def test_best_per_query_and_short_row_skipped(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(
                Path(td) / "blast.tsv",
                "q1\ts1\t90.0\t100\t95\t1e-10\t150\tsp1\n"
                "q1\ts2\t95.0\t100\t95\t1e-20\t200\tsp2\n"
                # Second hit for q1, worse bitscore than the current
                # best — must not overwrite it.
                "q1\ts4\t70.0\t100\t95\t1e-2\t10\tsp4\n"
                "\n"
                "q2\ts3\t80.0\t100\t95\t1e-5\t50\tsp3\n"
                "malformed\trow\n",
            )
            hits = collate.top_blast_hits(p)
            self.assertEqual(len(hits), 2)
            best_q1 = next(h for h in hits if h["qaccver"] == "q1")
            self.assertEqual(best_q1["saccver"], "s2")

    def test_homology_section_empty(self):
        self.assertIsNone(collate.homology_section(None))

    def test_homology_section_present(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(
                Path(td) / "blast.tsv",
                "q1\ts1\t90.0\t100\t95\t1e-10\t150\tsp1\n",
            )
            section = collate.homology_section(p)
            self.assertEqual(len(section["top_hits"]), 1)


class TestBarcodesSection(unittest.TestCase):
    def test_empty(self):
        self.assertIsNone(collate.barcodes_section(None))

    def test_present(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(
                Path(td) / "validation.tsv",
                "gene\tstatus\n COX1\tpass\nCOX2\tfail\n",
            )
            section = collate.barcodes_section(p)
            self.assertEqual(section["n_passed"], 1)

    def test_n_loci_passed(self):
        with tempfile.TemporaryDirectory() as td:
            p = write(
                Path(td) / "validation.tsv",
                "gene\tstatus\nCOX1\tpass\nCOX2\tpass\nCOX3\tnot_found\n",
            )
            self.assertEqual(collate.n_loci_passed(p), 2)
            self.assertEqual(collate.n_loci_passed(None), 0)


class TestWithCanonicalNames(unittest.TestCase):
    def _gene_sets_path(self, td):
        return write_json(Path(td) / "gene_sets.json", GENE_SETS)

    def test_no_crosscheck_returned_unchanged(self):
        summary = {"assembly_target": "animal_mt"}
        self.assertEqual(
            collate.with_canonical_names(summary, "/x"), summary)

    def test_no_gene_sets_path_returned_unchanged(self):
        summary = {
            "assembly_target": "animal_mt",
            "cds_crosscheck": {"agreed": ["cob_1"]},
        }
        self.assertEqual(
            collate.with_canonical_names(summary, None), summary)

    def test_transforms_agreed_and_miniprot_only(self):
        with tempfile.TemporaryDirectory() as td:
            gs = self._gene_sets_path(td)
            summary = {
                "assembly_target": "animal_mt",
                "cds_crosscheck": {
                    "agreed": ["cob_1"],
                    "miniprot_only": ["ATP8"],
                    "annotator_only": [
                        {"gene": "nad2_1", "reason": "off_panel"}],
                },
            }
            out = collate.with_canonical_names(summary, gs)
            xc = out["cds_crosscheck"]
            self.assertEqual(xc["agreed"]["raw"], ["cob_1"])
            self.assertEqual(xc["agreed"]["canonical"], ["CYTB"])
            self.assertEqual(xc["miniprot_only"]["canonical"], ["ATP8"])
            self.assertEqual(xc["annotator_only"][0]["canonical"], "NAD2")
            self.assertEqual(xc["annotator_only"][0]["reason"], "off_panel")

    def test_empty_annotator_only_left_alone(self):
        with tempfile.TemporaryDirectory() as td:
            gs = self._gene_sets_path(td)
            summary = {
                "assembly_target": "animal_mt",
                "cds_crosscheck": {"agreed": [], "annotator_only": []},
            }
            out = collate.with_canonical_names(summary, gs)
            self.assertEqual(out["cds_crosscheck"]["annotator_only"], [])


class TestProvenanceSection(unittest.TestCase):
    def test_no_annotation_summary_or_manifest(self):
        section = collate.provenance_section(None, None, None, "")
        self.assertIsNone(section["pipeline_commit"])
        self.assertIsNone(section["reference_bundle_version"])
        self.assertEqual(section["tool_versions"], collate.TOOL_VERSIONS)

    def test_with_annotation_summary_and_manifest(self):
        annotation_summary = {"tool_versions": {"miniprot": "0.18"}}
        manifest = {"version": "v2026.09", "generated_at": "2026-09-01"}
        section = collate.provenance_section(
            annotation_summary, {"selected_table": 5}, manifest, "abc123")
        self.assertEqual(section["pipeline_commit"], "abc123")
        self.assertEqual(section["reference_bundle_version"], "v2026.09")
        self.assertEqual(section["tool_versions"]["miniprot"], "0.18")
        self.assertEqual(section["genetic_code"], {"selected_table": 5})


class Args:
    """Minimal stand-in for argparse.Namespace, defaulting every
    optional collate.py CLI flag to None."""

    def __init__(self, **kwargs):
        defaults = dict(
            meta_json=None, pipeline_commit="",
            status_json=None, coverage_json=None,
            recruit_stats=None, nanoplot_raw=None, nanoplot_clean=None,
            bin_metadata_json=None, assembly_info=None,
            target_fasta=None, blast_tsv=None, barcodes_fasta=None,
            secondaries_tsv=None, validation_tsv=None,
            annotation_gff=None, annotation_summary=None,
            genetic_code_json=None, organelle_map_svg=None,
            graph_png=None, plastid_isoforms=None,
            gene_sets=None, refs_manifest=None, schema=None,
            params_json=None, report_templates=None, report_static=None,
            out_metadata=None, out_report=None,
        )
        defaults.update(kwargs)
        for key, value in defaults.items():
            setattr(self, key, value)


class TestBuildMetadataAndBundle(unittest.TestCase):
    """End-to-end cases against build_metadata + assemble_bundle,
    covering task 42 §8's numbered cases directly."""

    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.dir = Path(self.td.name)
        self.addCleanup(self.td.cleanup)
        self.cwd = Path.cwd()

    def _args(self, **kwargs):
        meta = kwargs.pop("meta", {
            "sample_id": "S1", "assembly_target": "animal_mt",
            "sample_info": "", "sample_type": "", "sample_receipt_date": "",
            "storage_location": "",
        })
        meta_path = write_json(self.dir / "meta.json", meta)
        status_json = kwargs.pop(
            "status_json", write_json(
                self.dir / "sample_status.json", {"status": "ok"}))
        coverage_json = kwargs.pop(
            "coverage_json", write_json(
                self.dir / "coverage.json", {"estimated_cov": 50.0}))
        return Args(
            meta_json=meta_path, status_json=status_json,
            coverage_json=coverage_json, **kwargs)

    def test_case3_low_coverage_is_full_bundle(self):
        args = self._args(
            status_json=write_json(
                self.dir / "s.json", {"status": "low_coverage"}),
            bin_metadata_json=write_json(
                self.dir / "b.json", {"contigs_selected": ["c1"]}),
            target_fasta=write(self.dir / "t.fasta", ">c1\nACGT\n"),
            validation_tsv=write(
                self.dir / "v.tsv", "gene\tstatus\nCOX1\tpass\n"),
        )
        metadata = collate.build_metadata(args)
        self.assertEqual(metadata["sample_status"], "low_coverage")
        self.assertEqual(metadata["bundle"], "full")

    def test_case4_no_barcode_is_full_bundle(self):
        args = self._args(
            bin_metadata_json=write_json(
                self.dir / "b.json", {"contigs_selected": ["c1"]}),
            target_fasta=write(self.dir / "t.fasta", ">c1\nACGT\n"),
            validation_tsv=write(
                self.dir / "v.tsv", "gene\tstatus\nCOX1\tnot_found\n"),
        )
        metadata = collate.build_metadata(args)
        self.assertEqual(metadata["sample_status"], "no_barcode")
        self.assertEqual(metadata["bundle"], "full")

    def test_case5_withheld_substitution_no_assembly(self):
        iso_dir = self.dir / "plastid_isoforms"
        iso_dir.mkdir()
        write(iso_dir / "path1.fasta", ">path1\nACGT\n")
        args = self._args(
            meta={
                "sample_id": "S1", "assembly_target": "plant_pt",
                "sample_info": "", "sample_type": "",
                "sample_receipt_date": "", "storage_location": "",
            },
            bin_metadata_json=write_json(
                self.dir / "b.json", {
                    "contigs_selected": [],
                    "plastid_canonicalisation": {
                        "branch": "canonical",
                        "substitution_applied": False,
                        "substitution_withheld_reason": "no_c3_selection",
                    },
                }),
            target_fasta=write(self.dir / "t.fasta", ""),
            plastid_isoforms=iso_dir,
        )
        metadata = collate.build_metadata(args)
        self.assertEqual(metadata["sample_status"], "no_assembly")
        self.assertEqual(metadata["bundle"], "minimal")

        import os
        os.chdir(self.dir)
        try:
            collate.assemble_bundle(args, metadata)
        finally:
            os.chdir(self.cwd)
        self.assertTrue((self.dir / "diagnostics" / "plastid_isoforms"
                         / "path1.fasta").is_file())
        self.assertFalse((self.dir / "organelle_assembly.fasta").exists())

    def test_case6_minimal_bundle_no_placeholders(self):
        args = self._args(
            status_json=write_json(self.dir / "s.json", {"status": "fail"}),
        )
        metadata = collate.build_metadata(args)
        self.assertEqual(metadata["bundle"], "minimal")

        import os
        os.chdir(self.dir)
        try:
            collate.assemble_bundle(args, metadata)
        finally:
            os.chdir(self.cwd)
        self.assertFalse((self.dir / "organelle_assembly.fasta").exists())
        self.assertFalse((self.dir / "organelle_annotation.gff").exists())
        self.assertFalse((self.dir / "barcodes.fasta").exists())
        self.assertFalse((self.dir / "diagnostics").exists())

    def test_case8_genetic_code_provenance_roundtrip(self):
        genetic_code = {
            "$schema": "wf5/genetic-code-selection/v1",
            "selected_table": 5,
            "candidates": [
                {"table": 2, "total_score": 10.0, "n_genes": 3},
                {"table": 5, "total_score": 20.0, "n_genes": 5},
            ],
        }
        args = self._args(
            bin_metadata_json=write_json(
                self.dir / "b.json", {"contigs_selected": ["c1"]}),
            target_fasta=write(self.dir / "t.fasta", ">c1\nACGT\n"),
            validation_tsv=write(
                self.dir / "v.tsv", "gene\tstatus\nCOX1\tpass\n"),
            genetic_code_json=write_json(
                self.dir / "gc.json", genetic_code),
        )
        metadata = collate.build_metadata(args)
        self.assertEqual(
            metadata["provenance"]["genetic_code"], genetic_code)

    def test_case9_missing_malformed_zero_byte_inputs(self):
        args = self._args(
            status_json=write(self.dir / "s.json", "{bad json"),
            bin_metadata_json=Path("/no/such/bin_metadata.json"),
            target_fasta=write(self.dir / "t.fasta", ""),
            annotation_summary=write(self.dir / "a.json", ""),
        )
        metadata = collate.build_metadata(args)
        # A malformed gate status.json parses to {} -> status None ->
        # unrecognised -> fails closed to "fail", never crashes.
        self.assertEqual(metadata["sample_status"], "fail")
        self.assertEqual(metadata["$schema"], collate.SCHEMA)

        import jsonschema
        jsonschema.validate(metadata, SCHEMA)

    def test_full_bundle_writes_all_five_files(self):
        args = self._args(
            bin_metadata_json=write_json(
                self.dir / "b.json", {"contigs_selected": ["c1"]}),
            target_fasta=write(self.dir / "t.fasta", ">c1\nACGT\n"),
            annotation_gff=write(self.dir / "a.gff", "##gff-version 3\n"),
            barcodes_fasta=write(self.dir / "bc.fasta", ">COX1\nACGT\n"),
            validation_tsv=write(
                self.dir / "v.tsv", "gene\tstatus\nCOX1\tpass\n"),
            secondaries_tsv=write(self.dir / "sec.tsv", "contig_id\n"),
            organelle_map_svg=write(self.dir / "m.svg", "<svg/>"),
            graph_png=write(self.dir / "g.png", "not-really-a-png"),
            annotation_summary=write_json(
                self.dir / "as.json", {
                    "tool_versions": {"miniprot": "0.18"},
                }),
        )
        metadata = collate.build_metadata(args)
        self.assertEqual(metadata["bundle"], "full")
        self.assertEqual(
            metadata["provenance"]["tool_versions"]["miniprot"], "0.18")

        import os
        os.chdir(self.dir)
        try:
            collate.assemble_bundle(args, metadata)
        finally:
            os.chdir(self.cwd)
        self.assertTrue((self.dir / "organelle_assembly.fasta").is_file())
        self.assertTrue((self.dir / "organelle_annotation.gff").is_file())
        self.assertTrue((self.dir / "barcodes.fasta").is_file())
        self.assertTrue((self.dir / "diagnostics" / "secondaries.tsv")
                        .is_file())
        self.assertTrue((self.dir / "diagnostics" / "graph.png").is_file())

        jsonschema_mod = __import__("jsonschema")
        jsonschema_mod.validate(metadata, SCHEMA)

    def test_assemble_bundle_full_with_no_source_files(self):
        """Defensive branch: `bundle == "full"` but the source paths
        for one or more of the three bundle files were never supplied
        — should not raise, and should simply skip that file."""
        args = self._args()
        import os
        os.chdir(self.dir)
        try:
            collate.assemble_bundle(args, {"bundle": "full"})
        finally:
            os.chdir(self.cwd)
        self.assertFalse((self.dir / "organelle_assembly.fasta").exists())
        self.assertFalse((self.dir / "organelle_annotation.gff").exists())
        self.assertFalse((self.dir / "barcodes.fasta").exists())

    def test_schema_valid_for_every_sample_status(self):
        import jsonschema

        cases = [
            ("fail", {"status_json": {"status": "fail"}}),
            ("no_assembly", {
                "bin_metadata_json": {"contigs_selected": []},
                "target_fasta": "",
            }),
            ("no_barcode", {
                "bin_metadata_json": {"contigs_selected": ["c1"]},
                "target_fasta": ">c1\nACGT\n",
                "validation_tsv": "gene\tstatus\nCOX1\tnot_found\n",
            }),
            ("low_coverage", {
                "status_json": {"status": "low_coverage"},
                "bin_metadata_json": {"contigs_selected": ["c1"]},
                "target_fasta": ">c1\nACGT\n",
                "validation_tsv": "gene\tstatus\nCOX1\tpass\n",
            }),
            ("ok", {
                "bin_metadata_json": {"contigs_selected": ["c1"]},
                "target_fasta": ">c1\nACGT\n",
                "validation_tsv": "gene\tstatus\nCOX1\tpass\n",
            }),
        ]
        for expected_status, overrides in cases:
            with self.subTest(status=expected_status):
                kwargs = {}
                if "status_json" in overrides:
                    kwargs["status_json"] = write_json(
                        self.dir / f"{expected_status}.status.json",
                        overrides["status_json"])
                if "bin_metadata_json" in overrides:
                    kwargs["bin_metadata_json"] = write_json(
                        self.dir / f"{expected_status}.bin.json",
                        overrides["bin_metadata_json"])
                if "target_fasta" in overrides:
                    kwargs["target_fasta"] = write(
                        self.dir / f"{expected_status}.fasta",
                        overrides["target_fasta"])
                if "validation_tsv" in overrides:
                    kwargs["validation_tsv"] = write(
                        self.dir / f"{expected_status}.tsv",
                        overrides["validation_tsv"])
                args = self._args(**kwargs)
                metadata = collate.build_metadata(args)
                self.assertEqual(metadata["sample_status"], expected_status)
                jsonschema.validate(metadata, SCHEMA)


class TestMain(unittest.TestCase):
    """CLI entry point — argument parsing, schema validation branches,
    and the always-exit-0 contract (task 42 §7)."""

    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.dir = Path(self.td.name)
        self.addCleanup(self.td.cleanup)
        self.cwd = Path.cwd()
        import os
        os.chdir(self.dir)
        self.addCleanup(os.chdir, self.cwd)

    def _run(self, extra_args):
        meta = write_json(self.dir / "meta.json", {
            "sample_id": "S1", "assembly_target": "animal_mt",
            "sample_info": "", "sample_type": "",
            "sample_receipt_date": "", "storage_location": "",
        })
        status = write_json(
            self.dir / "sample_status.json", {"status": "fail"})
        coverage = write_json(
            self.dir / "coverage.json", {"estimated_cov": 1.0})
        argv = [
            "collate.py",
            "--meta-json", str(meta),
            "--status-json", str(status),
            "--coverage-json", str(coverage),
        ] + extra_args
        old_argv = sys.argv
        sys.argv = argv
        try:
            return collate.main()
        finally:
            sys.argv = old_argv

    def test_exit_zero_and_outputs_written(self):
        rc = self._run([])
        self.assertEqual(rc, 0)
        self.assertTrue((self.dir / "metadata.json").is_file())
        self.assertTrue((self.dir / "report.html").is_file())

    def test_schema_flag_valid(self):
        schema_path = REPO_ROOT / "assets" / "sample_metadata.schema.json"
        rc = self._run(["--schema", str(schema_path)])
        self.assertEqual(rc, 0)

    def test_schema_flag_missing_file_skipped(self):
        rc = self._run(["--schema", "/no/such/schema.json"])
        self.assertEqual(rc, 0)

    def test_schema_flag_failing_validation_warns_not_raises(self):
        bad_schema = write_json(
            self.dir / "bad_schema.json",
            {"type": "object", "required": ["not_a_real_field"]})
        rc = self._run(["--schema", str(bad_schema)])
        self.assertEqual(rc, 0)


if __name__ == "__main__":
    unittest.main()
