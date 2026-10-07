"""Unit tests for bin/contig_bucket.py — task 52 §5.1.

The single shared classification -> bucket -> colour mapping consumed
by `bin/report/report.py` (contig table, coverage chart) and
`bin/annotate_graph_svg.py` (graph node colouring).

Cases:
  1. contig_bucket(): every classification bin_target.py actually
     emits, plus None and an unrecognised string.
  2. BUCKET_COLOURS keys == BUCKET_ORDER; MIXED_COLOUR is distinct
     from every bucket colour and from annotate_graph_svg's
     DISPLAY_COLOUR.
  3. bucket_legend(): shape and order.
  4. load_bin_metadata(): valid file, missing path, None, empty file,
     malformed JSON.
"""

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN_DIR))

import contig_bucket as cb  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "annotate_graph_svg", BIN_DIR / "annotate_graph_svg.py")
ags = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ags)


class ContigBucketTests(unittest.TestCase):

    def test_target_candidate(self):
        self.assertEqual(cb.contig_bucket("target_candidate"), "target")

    def test_secondary_target(self):
        self.assertEqual(cb.contig_bucket("secondary_target"), "secondary")

    def test_off_target(self):
        self.assertEqual(cb.contig_bucket("off_target"), "off-target")

    def test_sibling_organelle(self):
        self.assertEqual(cb.contig_bucket("sibling_organelle"), "off-target")

    def test_none(self):
        self.assertEqual(cb.contig_bucket(None), "off-target")

    def test_unrecognised_string(self):
        self.assertEqual(cb.contig_bucket("something_new"), "off-target")


class PaletteTests(unittest.TestCase):

    def test_bucket_colours_keys_match_order(self):
        self.assertEqual(set(cb.BUCKET_COLOURS), set(cb.BUCKET_ORDER))

    def test_mixed_colour_not_a_bucket_colour(self):
        self.assertNotIn(cb.MIXED_COLOUR, cb.BUCKET_COLOURS.values())

    def test_mixed_colour_not_display_colour(self):
        self.assertNotEqual(cb.MIXED_COLOUR, ags.DISPLAY_COLOUR)


class BucketLegendTests(unittest.TestCase):

    def test_shape_and_order(self):
        legend = cb.bucket_legend()
        self.assertEqual(
            [row["key"] for row in legend], cb.BUCKET_ORDER)
        for row in legend:
            self.assertEqual(row["name"], cb.BUCKET_LABELS[row["key"]])
            self.assertEqual(row["color"], cb.BUCKET_COLOURS[row["key"]])


class LoadBinMetadataTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_valid_file(self):
        path = self.tmp / "bin_metadata.json"
        path.write_text(json.dumps({"contigs": [{"contig_id": "c1"}]}))
        self.assertEqual(
            cb.load_bin_metadata(path), {"contigs": [{"contig_id": "c1"}]})

    def test_missing_path(self):
        self.assertEqual(
            cb.load_bin_metadata(self.tmp / "absent.json"), {})

    def test_none(self):
        self.assertEqual(cb.load_bin_metadata(None), {})

    def test_empty_file(self):
        path = self.tmp / "empty.json"
        path.write_text("")
        self.assertEqual(cb.load_bin_metadata(path), {})

    def test_malformed_json(self):
        path = self.tmp / "bad.json"
        path.write_text("{not json")
        self.assertEqual(cb.load_bin_metadata(path), {})


if __name__ == "__main__":
    unittest.main()
