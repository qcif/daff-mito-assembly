"""Unit tests for bin/annotation_scores.py — task 40 §6.

Cases (numbered to match task 40 §6):
  1. Best-hit selection when a query has multiple HSPs and multiple
     subjects — highest bitscore wins; a bitscore tie breaks on
     identity.
  2. A query with no blastp hit at all — explicit null triplet, the
     feature still appears in the GFF and in cds_scores (not zero,
     not dropped).
  3. Features from both sources (miniprot mRNA, rescued mRNA) scored
     identically by the same code path.
  4. Malformed/truncated blastp output — skipped with a warning, the
     rest of the file still processed.
  5. The no-annotator targets (plant_pt/plant_mt) — miniprot-only
     features score normally with no mitos features present at all.
  6. Non-CDS (tRNA/rRNA) and non-mRNA lines pass through untouched.
  7. Missing blast.tsv file entirely — every feature gets an explicit
     null triplet, not an error.
  8. $schema is bumped to v3 on the written summary.

Plus: blank lines in blast.tsv skipped, an mRNA row with no ID=
attribute left untouched, and main() wired end-to-end for 100%
branch coverage (rule 14).
"""

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN_DIR))

_spec = importlib.util.spec_from_file_location(
    "annotation_scores", BIN_DIR / "annotation_scores.py")
ascores = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ascores)


MRNA_MP1 = (
    "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\tID=MP1;Target=panel_ATP6 1 30"
)
CDS_MP1 = "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP1"
MRNA_R1 = (
    "ctg1\tmitos\tmRNA\t100\t150\t.\t+\t.\tID=r1;Name=ATP8"
)
CDS_R1 = "ctg1\tmitos\tCDS\t100\t150\t.\t+\t0\tParent=r1"
TRNA_LINE = (
    "ctg1\tmitos\ttRNA\t200\t270\t.\t+\t.\tID=trn1;Name=trnF(gaa)"
)


class BaseCase(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _write(self, name, text):
        path = self.tmp / name
        path.write_text(text)
        return path

    def _run(self, gff_lines, blast_lines, summary=None):
        gff_path = self._write("merged.gff", "\n".join(gff_lines) + "\n")
        blast_path = self._write(
            "blast.tsv",
            "\n".join(blast_lines) + ("\n" if blast_lines else ""))
        summary_path = self._write("annotation_summary.json", json.dumps(
            summary if summary is not None else {"$schema": "old"}))
        out_gff = self.tmp / "out.gff"
        out_summary = self.tmp / "out_summary.json"
        ascores.run(
            gff_path, summary_path, blast_path, out_gff, out_summary)
        return (
            out_gff.read_text(),
            json.loads(out_summary.read_text()),
        )


class AnnotationScoresTests(BaseCase):

    def test_best_hit_by_bitscore_multiple_hsps_and_subjects(self):
        gff_text, summary = self._run(
            [MRNA_MP1, CDS_MP1],
            [
                "MP1\tref_a\t50.0\t60.0\t1e-5\t100.0",
                "MP1\tref_b\t95.0\t99.0\t1e-40\t250.0",
                "MP1\tref_c\t70.0\t80.0\t1e-10\t150.0",
            ],
        )
        row = summary["cds_scores"][0]
        self.assertEqual(row["bitscore"], 250.0)
        self.assertEqual(row["pident"], 95.0)
        self.assertEqual(row["sseqid"], "ref_b")
        self.assertIn("pident=95.0;qcovhsp=99.0;bitscore=250.0", gff_text)

    def test_bitscore_tie_broken_on_identity(self):
        _gff, summary = self._run(
            [MRNA_MP1, CDS_MP1],
            [
                "MP1\tref_a\t80.0\t90.0\t1e-10\t200.0",
                "MP1\tref_b\t92.0\t90.0\t1e-10\t200.0",
            ],
        )
        self.assertEqual(summary["cds_scores"][0]["sseqid"], "ref_b")

    def test_no_hit_at_all_null_triplet_feature_kept(self):
        gff_text, summary = self._run([MRNA_MP1, CDS_MP1], [])
        row = summary["cds_scores"][0]
        self.assertIsNone(row["pident"])
        self.assertIsNone(row["qcovhsp"])
        self.assertIsNone(row["bitscore"])
        self.assertIn("ID=MP1", gff_text)
        self.assertIn("pident=null;qcovhsp=null;bitscore=null", gff_text)

    def test_miniprot_and_rescued_scored_identically(self):
        _gff, summary = self._run(
            [MRNA_MP1, CDS_MP1, MRNA_R1, CDS_R1],
            [
                "MP1\tref_a\t80.0\t90.0\t1e-10\t200.0",
                "r1\tref_b\t70.0\t85.0\t1e-8\t150.0",
            ],
        )
        by_id = {row["id"]: row for row in summary["cds_scores"]}
        self.assertEqual(by_id["MP1"]["bitscore"], 200.0)
        self.assertEqual(by_id["r1"]["bitscore"], 150.0)
        self.assertEqual(by_id["r1"]["gene"], "ATP8")

    def test_malformed_blast_lines_skipped(self):
        gff_text, summary = self._run(
            [MRNA_MP1, CDS_MP1],
            [
                "too\tfew\tcolumns",
                "MP1\tref_a\tNOT_A_NUMBER\t90.0\t1e-10\t200.0",
                "MP1\tref_b\t80.0\t90.0\t1e-10\t200.0",
            ],
        )
        row = summary["cds_scores"][0]
        self.assertEqual(row["bitscore"], 200.0)
        self.assertEqual(row["sseqid"], "ref_b")

    def test_no_annotator_target_miniprot_only(self):
        gff_text, summary = self._run(
            [MRNA_MP1, CDS_MP1],
            ["MP1\tref_a\t80.0\t90.0\t1e-10\t200.0"],
        )
        self.assertEqual(len(summary["cds_scores"]), 1)
        self.assertEqual(summary["cds_scores"][0]["source"], "miniprot")

    def test_non_cds_lines_pass_through_untouched(self):
        gff_text, _summary = self._run(
            [MRNA_MP1, CDS_MP1, TRNA_LINE],
            ["MP1\tref_a\t80.0\t90.0\t1e-10\t200.0"],
        )
        self.assertIn(TRNA_LINE, gff_text)

    def test_missing_blast_tsv_file_all_null(self):
        blast_path = self.tmp / "does_not_exist.tsv"
        gff_path = self._write(
            "merged.gff", "\n".join([MRNA_MP1, CDS_MP1]) + "\n")
        summary_path = self._write(
            "annotation_summary.json", json.dumps({}))
        out_gff = self.tmp / "out.gff"
        out_summary = self.tmp / "out_summary.json"
        ascores.run(
            gff_path, summary_path, blast_path, out_gff, out_summary)
        summary = json.loads(out_summary.read_text())
        self.assertIsNone(summary["cds_scores"][0]["bitscore"])

    def test_schema_bumped_to_v3(self):
        _gff, summary = self._run([MRNA_MP1, CDS_MP1], [])
        self.assertEqual(summary["$schema"], "wf5/annotation-summary/v3")

    def test_blank_blast_lines_skipped(self):
        _gff, summary = self._run(
            [MRNA_MP1, CDS_MP1],
            ["", "MP1\tref_a\t80.0\t90.0\t1e-10\t200.0", ""],
        )
        self.assertEqual(summary["cds_scores"][0]["bitscore"], 200.0)

    def test_mrna_without_id_left_untouched(self):
        line = "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\tTarget=panel_ATP6 1 30"
        result = ascores.rescored_line(line, {})
        self.assertEqual(result, line)

    def test_main_happy_path(self):
        gff_path = self._write(
            "merged.gff", "\n".join([MRNA_MP1, CDS_MP1]) + "\n")
        summary_path = self._write(
            "annotation_summary.json", json.dumps({}))
        blast_path = self._write(
            "blast.tsv", "MP1\tref_a\t80.0\t90.0\t1e-10\t200.0\n")
        out_gff = self.tmp / "out.gff"
        out_summary = self.tmp / "out_summary.json"

        old_argv = sys.argv
        sys.argv = [
            "annotation_scores.py",
            "--merged-gff", str(gff_path),
            "--annotation-summary", str(summary_path),
            "--blast-tsv", str(blast_path),
            "--out-gff", str(out_gff),
            "--out-summary", str(out_summary),
        ]
        try:
            rc = ascores.main()
        finally:
            sys.argv = old_argv
        self.assertEqual(rc, 0)
        self.assertTrue(out_summary.exists())


if __name__ == "__main__":
    unittest.main()
