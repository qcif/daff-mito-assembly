"""Unit tests for bin/annotation_gff.py — task 40 §2 (shared parser).

Cases:
  1. Miniprot-provenance mRNA (Target=, no Name=) — gene from Target.
  2. Rescued (mitos)-provenance mRNA (Name=) — gene from Name.
  3. Multi-exon feature — CDS intervals accumulate onto its parent.
  4. mRNA with no CDS children — dropped from the result.
  5. CDS row with unknown Parent — ignored, no KeyError.
  6. mRNA row missing ID= — skipped.
  7. Non-9-column / comment / blank lines — skipped, not fatal.
  8. Non-numeric start/end on mRNA or CDS — skipped, not fatal.
  9. parse_attrs: entries without '=' are ignored.
"""

import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[2] / "bin" / "annotation_gff.py"
)
_spec = importlib.util.spec_from_file_location(
    "annotation_gff", MODULE_PATH)
agff = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(agff)


class ParseAnnotationGffTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _write(self, lines):
        path = self.tmp / "merged.gff"
        path.write_text("\n".join(lines) + "\n")
        return path

    def test_miniprot_gene_from_target(self):
        path = self._write([
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\t"
            "ID=MP1;Target=panel_ATP6 1 30",
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP1",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]["gene"], "ATP6")
        self.assertEqual(features[0]["source"], "miniprot")
        self.assertEqual(features[0]["cds"], [(1, 90)])

    def test_rescued_gene_from_name(self):
        path = self._write([
            "ctg1\tmitos\tmRNA\t1\t50\t.\t+\t.\t"
            "ID=mitos_rescue_ATP8_1;Name=ATP8;OriginalName=atp8_0",
            "ctg1\tmitos\tCDS\t1\t50\t.\t+\t0\t"
            "ID=mitos_rescue_ATP8_1_cds1;Parent=mitos_rescue_ATP8_1",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]["gene"], "ATP8")
        self.assertEqual(features[0]["source"], "mitos")

    def test_multi_exon_accumulates_cds(self):
        path = self._write([
            "ctg1\tmitos\tmRNA\t1\t100\t.\t+\t.\tID=r1;Name=NAD1",
            "ctg1\tmitos\tCDS\t1\t40\t.\t+\t0\tParent=r1",
            "ctg1\tmitos\tCDS\t60\t100\t.\t+\t1\tParent=r1",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(features[0]["cds"], [(1, 40), (60, 100)])

    def test_mrna_without_cds_dropped(self):
        path = self._write([
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\t"
            "ID=MP1;Target=panel_ATP6 1 30",
        ])
        self.assertEqual(agff.parse_annotation_gff(path), [])

    def test_cds_with_unknown_parent_ignored(self):
        path = self._write([
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=ghost",
        ])
        self.assertEqual(agff.parse_annotation_gff(path), [])

    def test_mrna_missing_id_skipped(self):
        path = self._write([
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\tTarget=panel_ATP6 1 30",
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP1",
        ])
        self.assertEqual(agff.parse_annotation_gff(path), [])

    def test_malformed_lines_skipped(self):
        path = self._write([
            "##gff-version 3",
            "",
            "ctg1\tminiprot\tmRNA\ttoo\tfew\tcolumns",
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\t"
            "ID=MP1;Target=panel_ATP6 1 30",
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP1",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(len(features), 1)

    def test_non_numeric_coordinates_skipped(self):
        path = self._write([
            "ctg1\tminiprot\tmRNA\tNA\tNA\t.\t+\t.\t"
            "ID=MP1;Target=panel_ATP6 1 30",
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\t"
            "ID=MP2;Target=panel_ATP6 1 30",
            "ctg1\tminiprot\tCDS\tNA\tNA\t.\t+\t0\tParent=MP2",
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP2",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]["id"], "MP2")
        self.assertEqual(features[0]["cds"], [(1, 90)])

    def test_parse_attrs_ignores_entries_without_equals(self):
        attrs = agff.parse_attrs("ID=x;standalone;Name=y")
        self.assertEqual(attrs, {"ID": "x", "Name": "y"})

    def test_gene_symbol_no_target_no_name(self):
        self.assertEqual(agff.gene_symbol({}), "")


if __name__ == "__main__":
    unittest.main()
