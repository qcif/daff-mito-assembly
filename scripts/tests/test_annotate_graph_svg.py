"""Unit tests for bin/annotate_graph_svg.py — task 47 (stage 10 C12).

Fixtures are small hand-written GFA/assembly_info.txt/SVG fragments
shaped like the real BandageNG/Flye output captured in task 47 §2, not
captured tool output itself (CONSTITUTION rule 19).

Cases:
  1. parse_gfa_segments: `S` records in order; non-`S` lines and short
     lines skipped.
  2. allocate_sentinels / write_sentinel_csv: deterministic, unique,
     never collides with white/black, header row present.
  3. run_allocate: GFA in, CSV out, end to end.
  4. parse_graph_path_column: missing file; signed entries fold;
     `*` skipped; comment lines skipped; malformed (<8-field) lines
     skipped; one segment traversed by two contigs; multiple segments.
  5. read_sentinel_csv: header skipped, malformed rows skipped.
  6. annotate_svg: sentinel match wraps the path with data-node/
     data-contigs/<title> and recolours the fill; a segment with no
     traversing contig gets an empty data-contigs and a bare title; a
     sentinel absent from the SVG warns; a path fill matching no
     sentinel warns unless it is the white background or black stroke;
     default no-op `warn`.
  7. strip_xml_prolog: XML declaration stripped, DOCTYPE stripped,
     both stripped together, absent prolog left alone.
  8. run_annotate: end-to-end file read/write, including prolog
     stripping on real BandageNG-shaped input.
  9. main(): `allocate` subcommand; `annotate` subcommand normal path;
     `annotate` exception path with the input SVG present (untouched
     copy) and absent (empty output) — always exits 0.
  10. segment_bucket (task 52): single-bucket cases, two/three-way
      mixed, an unclassified contig ignored, all-unclassified -> None,
      empty contig_ids -> None.
  11. annotate_svg bucket colouring (task 52): fill/data-bucket/title
      suffix per bucket, mixed, and unknown (no bin_metadata); unknown
      node's <title> stays byte-identical to task 47's.
  12. run_annotate/main --bin-metadata (task 52): omitted, missing
      file, empty file, malformed JSON all degrade to the all-neutral
      output and exit 0.
"""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN_DIR))

import contig_bucket  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "annotate_graph_svg", BIN_DIR / "annotate_graph_svg.py")
ags = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ags)


def _write(tmp: Path, name: str, text: str) -> Path:
    path = tmp / name
    path.write_text(text)
    return path


GFA_3_SEGMENTS = (
    "H\tVN:Z:1.0\n"
    "S\tedge_1\t*\tLN:i:100\n"
    "L\tedge_1\t+\tedge_2\t+\t0M\n"
    "S\tedge_2\t*\tLN:i:200\n"
    "S\tedge_3\t*\tLN:i:300\n"
)

# Shaped like the real BandageNG Qt SVG export (verified against the
# pinned biocontainer, task 47 §2/§4): fill lives on the <g> wrapper,
# not the <path> it contains, and background/stroke groups never carry
# a fill attribute at all.
RAW_SVG_3_SEGMENTS = (
    '<svg xmlns="http://www.w3.org/2000/svg">'
    '<desc>Generated with Qt</desc>'
    '<g fill="#ffffff"><path d="M 0 0 L 10 10 Z"/></g>'
    '<g fill="#000001" stroke="none">'
    '<path vector-effect="none" d="M 1 1 L 2 2 Z"/></g>'
    '<g fill="#000002" stroke="none">'
    '<path vector-effect="none" d="M 3 3 L 4 4 Z"/></g>'
    '<g fill="#000003" stroke="none">'
    '<path vector-effect="none" d="M 5 5 L 6 6 Z"/></g>'
    '<g stroke="#000000"><path d="M 7 7 L 8 8 Z"/></g>'
    '</svg>'
)

ASSEMBLY_INFO = (
    "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
    "contig_1\t500\t20\tY\tN\t1\t*\t-2,1,2\n"
    "contig_2\t300\t15\tN\tN\t1\t*\t*,3,*\n"
)


class ParseGfaSegmentsTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_s_records_in_order(self):
        gfa = _write(self.tmp, "a.gfa", GFA_3_SEGMENTS)
        self.assertEqual(
            ags.parse_gfa_segments(gfa), ["edge_1", "edge_2", "edge_3"])

    def test_short_lines_skipped(self):
        gfa = _write(self.tmp, "a.gfa", "S\n")
        self.assertEqual(ags.parse_gfa_segments(gfa), [])

    def test_zero_segments(self):
        gfa = _write(self.tmp, "a.gfa", "H\tVN:Z:1.0\n")
        self.assertEqual(ags.parse_gfa_segments(gfa), [])


class AllocateSentinelsTests(unittest.TestCase):

    def test_deterministic(self):
        segments = ["edge_1", "edge_2", "edge_3"]
        self.assertEqual(
            ags.allocate_sentinels(segments), ags.allocate_sentinels(segments))

    def test_unique_and_no_collision(self):
        sentinels = ags.allocate_sentinels(["edge_1", "edge_2", "edge_3"])
        colours = [c for _, c in sentinels]
        self.assertEqual(len(colours), len(set(colours)))
        self.assertNotIn("#ffffff", colours)
        self.assertNotIn("#000000", colours)

    def test_one_segment(self):
        self.assertEqual(
            ags.allocate_sentinels(["edge_1"]), [("edge_1", "#000001")])


class WriteSentinelCsvTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_header_row_present(self):
        out = self.tmp / "sentinels.csv"
        ags.write_sentinel_csv([("edge_1", "#000001")], out)
        lines = out.read_text().splitlines()
        self.assertEqual(lines[0], "name,color")
        self.assertEqual(lines[1], "edge_1,#000001")


class RunAllocateTests(unittest.TestCase):

    def test_end_to_end(self):
        tmp = Path(tempfile.mkdtemp())
        gfa = _write(tmp, "a.gfa", GFA_3_SEGMENTS)
        out = tmp / "sentinels.csv"
        ags.run_allocate(gfa, out)
        lines = out.read_text().splitlines()
        self.assertEqual(lines[0], "name,color")
        self.assertEqual(len(lines), 4)


class ParseGraphPathColumnTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_missing_file(self):
        self.assertEqual(ags.parse_graph_path_column(None), {})
        self.assertEqual(
            ags.parse_graph_path_column(self.tmp / "absent.txt"), {})

    def test_signed_entries_fold_and_star_skipped(self):
        info = _write(self.tmp, "assembly_info.txt", ASSEMBLY_INFO)
        traversal = ags.parse_graph_path_column(info)
        self.assertEqual(traversal["edge_2"], {"contig_1"})
        self.assertEqual(traversal["edge_1"], {"contig_1"})
        self.assertEqual(traversal["edge_3"], {"contig_2"})

    def test_two_contigs_traverse_one_segment(self):
        info = _write(self.tmp, "assembly_info.txt", (
            "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\t"
            "graph_path\n"
            "contig_1\t500\t20\tY\tN\t1\t*\t1\n"
            "contig_2\t300\t15\tN\tN\t1\t*\t-1\n"
        ))
        self.assertEqual(
            ags.parse_graph_path_column(info)["edge_1"],
            {"contig_1", "contig_2"})

    def test_comment_and_malformed_lines_skipped(self):
        info = _write(self.tmp, "assembly_info.txt", (
            "#comment\n"
            "too\tfew\tcolumns\n"
            "contig_1\t500\t20\tY\tN\t1\t*\tnot_a_number\n"
        ))
        self.assertEqual(ags.parse_graph_path_column(info), {})


class ReadSentinelCsvTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_header_skipped_and_rows_parsed(self):
        csv = _write(self.tmp, "s.csv", "name,color\nedge_1,#000001\n")
        self.assertEqual(ags.read_sentinel_csv(csv), [("edge_1", "#000001")])

    def test_malformed_rows_skipped(self):
        csv = _write(self.tmp, "s.csv", "name,color\nedge_1,#000001,extra\n")
        self.assertEqual(ags.read_sentinel_csv(csv), [])


class AnnotateSvgTests(unittest.TestCase):

    def setUp(self):
        self.sentinels = [
            ("edge_1", "#000001"),
            ("edge_2", "#000002"),
            ("edge_3", "#000003"),
        ]
        self.traversal = {"edge_1": {"contig_1"}, "edge_2": {"contig_1"}}

    def test_matched_sentinels_wrapped_and_recoloured(self):
        warnings = []
        result = ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, self.sentinels, self.traversal,
            warn=warnings.append)
        self.assertIn('data-node="edge_1"', result)
        self.assertIn('data-contigs="contig_1"', result)
        self.assertIn('<title>edge_1 — contig_1</title>', result)
        self.assertIn(f'fill="{ags.DISPLAY_COLOUR}"', result)
        self.assertNotIn('fill="#000001"', result)
        self.assertNotIn('fill="#000002"', result)
        self.assertNotIn('fill="#000003"', result)

    def test_segment_with_no_traversing_contig(self):
        result = ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, self.sentinels, self.traversal)
        self.assertIn('data-node="edge_3"', result)
        self.assertIn('data-contigs=""', result)
        self.assertIn('<title>edge_3</title>', result)

    def test_sentinel_not_drawn_warns(self):
        warnings = []
        sentinels = self.sentinels + [("edge_4", "#000004")]
        ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, sentinels, self.traversal,
            warn=warnings.append)
        self.assertTrue(
            any("#000004" in w and "edge_4" in w for w in warnings))

    def test_unmatched_path_fill_warns(self):
        warnings = []
        svg = RAW_SVG_3_SEGMENTS.replace(
            'fill="#000003"', 'fill="#abcdef"')
        ags.annotate_svg(
            svg, self.sentinels, self.traversal, warn=warnings.append)
        self.assertTrue(any("#abcdef" in w for w in warnings))

    def test_background_and_stroke_colours_never_warn(self):
        warnings = []
        ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, self.sentinels, self.traversal,
            warn=warnings.append)
        self.assertFalse(any("#ffffff" in w for w in warnings))
        self.assertFalse(any("#000000" in w for w in warnings))

    def test_default_warn_is_noop(self):
        # No warn kwarg — the default no-op path must not raise, even
        # though this input has an unmatched sentinel.
        sentinels = self.sentinels + [("edge_4", "#000004")]
        result = ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, sentinels, self.traversal)
        self.assertIn('data-node="edge_1"', result)


class StripXmlPrologTests(unittest.TestCase):

    def test_xml_declaration_stripped(self):
        svg = '<?xml version="1.0" encoding="UTF-8"?>\n<svg></svg>'
        self.assertEqual(ags.strip_xml_prolog(svg), '<svg></svg>')

    def test_doctype_stripped(self):
        svg = '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN">\n<svg></svg>'
        self.assertEqual(ags.strip_xml_prolog(svg), '<svg></svg>')

    def test_both_stripped_together(self):
        svg = (
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN">\n'
            '<svg></svg>'
        )
        self.assertEqual(ags.strip_xml_prolog(svg), '<svg></svg>')

    def test_no_prolog_left_alone(self):
        self.assertEqual(ags.strip_xml_prolog('<svg></svg>'), '<svg></svg>')


class RunAnnotateTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_end_to_end(self):
        svg = _write(self.tmp, "raw.svg", RAW_SVG_3_SEGMENTS)
        csv = self.tmp / "sentinels.csv"
        ags.write_sentinel_csv(
            [("edge_1", "#000001"), ("edge_2", "#000002"),
             ("edge_3", "#000003")], csv)
        info = _write(self.tmp, "assembly_info.txt", ASSEMBLY_INFO)
        out = self.tmp / "graph.svg"
        ags.run_annotate(svg, csv, info, out, "SAMPLE01")
        self.assertIn('data-node="edge_1"', out.read_text())

    def test_xml_prolog_stripped_from_output(self):
        svg = _write(self.tmp, "raw.svg", (
            '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
            + RAW_SVG_3_SEGMENTS
        ))
        csv = self.tmp / "sentinels.csv"
        ags.write_sentinel_csv([("edge_1", "#000001")], csv)
        out = self.tmp / "graph.svg"
        ags.run_annotate(svg, csv, None, out, "SAMPLE01")
        result = out.read_text()
        self.assertNotIn('<?xml', result)
        self.assertTrue(result.startswith('<svg'))

    def test_zero_byte_svg_left_empty(self):
        svg = _write(self.tmp, "raw.svg", "")
        csv = self.tmp / "sentinels.csv"
        ags.write_sentinel_csv([("edge_1", "#000001")], csv)
        out = self.tmp / "graph.svg"
        ags.run_annotate(svg, csv, None, out, "SAMPLE01")
        self.assertEqual(out.read_text(), "")


class MainTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_allocate_subcommand(self):
        gfa = _write(self.tmp, "a.gfa", GFA_3_SEGMENTS)
        out = self.tmp / "sentinels.csv"
        argv = [
            "annotate_graph_svg.py", "allocate",
            "--gfa", str(gfa), "--out", str(out),
        ]
        with patch.object(sys, "argv", argv):
            self.assertEqual(ags.main(), 0)
        self.assertTrue(out.is_file())

    def test_annotate_subcommand_normal(self):
        svg = _write(self.tmp, "raw.svg", RAW_SVG_3_SEGMENTS)
        csv = self.tmp / "sentinels.csv"
        ags.write_sentinel_csv(
            [("edge_1", "#000001"), ("edge_2", "#000002"),
             ("edge_3", "#000003")], csv)
        info = _write(self.tmp, "assembly_info.txt", ASSEMBLY_INFO)
        out = self.tmp / "graph.svg"
        argv = [
            "annotate_graph_svg.py", "annotate",
            "--svg", str(svg), "--sentinels", str(csv),
            "--assembly-info", str(info), "--sample-id", "SAMPLE01",
            "--out", str(out),
        ]
        with patch.object(sys, "argv", argv):
            self.assertEqual(ags.main(), 0)
        self.assertIn('data-node="edge_1"', out.read_text())

    def test_annotate_subcommand_exception_svg_present_left_untouched(self):
        svg = _write(self.tmp, "raw.svg", RAW_SVG_3_SEGMENTS)
        missing_csv = self.tmp / "absent.csv"
        out = self.tmp / "graph.svg"
        argv = [
            "annotate_graph_svg.py", "annotate",
            "--svg", str(svg), "--sentinels", str(missing_csv),
            "--sample-id", "SAMPLE01", "--out", str(out),
        ]
        with patch.object(sys, "argv", argv):
            self.assertEqual(ags.main(), 0)
        self.assertEqual(out.read_text(), RAW_SVG_3_SEGMENTS)

    def test_annotate_subcommand_exception_svg_absent_writes_empty(self):
        missing_svg = self.tmp / "absent.svg"
        missing_csv = self.tmp / "absent.csv"
        out = self.tmp / "graph.svg"
        argv = [
            "annotate_graph_svg.py", "annotate",
            "--svg", str(missing_svg), "--sentinels", str(missing_csv),
            "--sample-id", "SAMPLE01", "--out", str(out),
        ]
        with patch.object(sys, "argv", argv):
            self.assertEqual(ags.main(), 0)
        self.assertEqual(out.read_text(), "")


class SegmentBucketTests(unittest.TestCase):

    def setUp(self):
        self.classifications = {
            "contig_1": "target_candidate",
            "contig_2": "secondary_target",
            "contig_3": "off_target",
        }

    def test_target_only(self):
        self.assertEqual(
            ags.segment_bucket({"contig_1"}, self.classifications),
            "target")

    def test_secondary_only(self):
        self.assertEqual(
            ags.segment_bucket({"contig_2"}, self.classifications),
            "secondary")

    def test_off_target_only(self):
        self.assertEqual(
            ags.segment_bucket({"contig_3"}, self.classifications),
            "off-target")

    def test_two_way_mixed(self):
        self.assertEqual(
            ags.segment_bucket(
                {"contig_1", "contig_3"}, self.classifications),
            "mixed")

    def test_three_way_mixed(self):
        self.assertEqual(
            ags.segment_bucket(
                {"contig_1", "contig_2", "contig_3"}, self.classifications),
            "mixed")

    def test_unclassified_contig_ignored(self):
        self.assertEqual(
            ags.segment_bucket(
                {"contig_1", "contig_x"}, self.classifications),
            "target")

    def test_all_unclassified_returns_none(self):
        self.assertIsNone(
            ags.segment_bucket({"contig_x", "contig_y"}, self.classifications))

    def test_empty_contig_ids_returns_none(self):
        self.assertIsNone(ags.segment_bucket(set(), self.classifications))


class AnnotateSvgBucketColouringTests(unittest.TestCase):
    """task 52 §3.4 -- bucket colouring layered onto annotate_svg()."""

    def setUp(self):
        self.sentinels = [
            ("edge_1", "#000001"),
            ("edge_2", "#000002"),
            ("edge_3", "#000003"),
        ]
        # edge_1: target only. edge_2: target + off-target (mixed).
        # edge_3: no traversing contig at all (unknown).
        self.traversal = {
            "edge_1": {"contig_1"},
            "edge_2": {"contig_1", "contig_2"},
        }
        self.classifications = {
            "contig_1": "target_candidate",
            "contig_2": "off_target",
        }

    def test_single_bucket_node_coloured_and_tagged(self):
        result = ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, self.sentinels, self.traversal,
            self.classifications)
        self.assertIn(
            f'fill="{contig_bucket.BUCKET_COLOURS["target"]}"', result)
        self.assertIn('data-bucket="target"', result)
        self.assertIn('<title>edge_1 — contig_1 (target)</title>', result)

    def test_mixed_node_gets_mixed_colour_and_ordered_tooltip(self):
        result = ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, self.sentinels, self.traversal,
            self.classifications)
        self.assertIn(f'fill="{contig_bucket.MIXED_COLOUR}"', result)
        self.assertIn('data-bucket="mixed"', result)
        # Fixed target/secondary/off-target order, not alphabetical or
        # set-iteration order.
        self.assertIn(
            '<title>edge_2 — contig_1,contig_2 '
            '(mixed: target, off-target)</title>',
            result)

    def test_unknown_node_matches_task_47_output_exactly(self):
        # No bin_metadata at all (contig_classifications omitted) ->
        # every node must render byte-identical to task 47's output.
        no_bucket_result = ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, self.sentinels, self.traversal)
        legacy_result = ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, self.sentinels, self.traversal, {})
        self.assertEqual(no_bucket_result, legacy_result)
        self.assertIn(f'fill="{ags.DISPLAY_COLOUR}"', legacy_result)
        self.assertNotIn('data-bucket=', legacy_result)
        self.assertIn('<title>edge_3</title>', legacy_result)

    def test_segment_with_no_traversing_contig_is_unknown(self):
        result = ags.annotate_svg(
            RAW_SVG_3_SEGMENTS, self.sentinels, self.traversal,
            self.classifications)
        self.assertIn(f'fill="{ags.DISPLAY_COLOUR}"', result)
        self.assertIn('<title>edge_3</title>', result)


BIN_METADATA_JSON = (
    '{"contigs": ['
    '{"contig_id": "contig_1", "classification": "target_candidate"}, '
    '{"contig_id": "contig_2", "classification": "off_target"}'
    ']}'
)


class RunAnnotateBinMetadataTests(unittest.TestCase):
    """task 52 §3.4/§5.2 -- --bin-metadata loading and its degrade
    path. A missing/empty/malformed file is not an exception: it is
    the "unknown" (all-neutral) row, i.e. today's output exactly."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.svg = _write(self.tmp, "raw.svg", RAW_SVG_3_SEGMENTS)
        self.csv = self.tmp / "sentinels.csv"
        ags.write_sentinel_csv(
            [("edge_1", "#000001"), ("edge_2", "#000002"),
             ("edge_3", "#000003")], self.csv)
        self.info = _write(self.tmp, "assembly_info.txt", (
            "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\t"
            "graph_path\n"
            "contig_1\t500\t20\tY\tN\t1\t*\t1\n"
            "contig_2\t300\t15\tN\tN\t1\t*\t2\n"
        ))

    def test_bin_metadata_omitted_is_all_neutral(self):
        out = self.tmp / "graph.svg"
        ags.run_annotate(
            self.svg, self.csv, self.info, out, "SAMPLE01")
        result = out.read_text()
        self.assertIn(f'fill="{ags.DISPLAY_COLOUR}"', result)
        self.assertNotIn('data-bucket=', result)

    def test_bin_metadata_present_colours_by_bucket(self):
        meta = _write(self.tmp, "bin_metadata.json", BIN_METADATA_JSON)
        out = self.tmp / "graph.svg"
        ags.run_annotate(
            self.svg, self.csv, self.info, out, "SAMPLE01", meta)
        result = out.read_text()
        self.assertIn('data-bucket="target"', result)

    def test_bin_metadata_missing_file_degrades(self):
        out = self.tmp / "graph.svg"
        ags.run_annotate(
            self.svg, self.csv, self.info, out, "SAMPLE01",
            self.tmp / "absent.json")
        result = out.read_text()
        self.assertNotIn('data-bucket=', result)

    def test_bin_metadata_empty_file_degrades(self):
        meta = _write(self.tmp, "bin_metadata.json", "")
        out = self.tmp / "graph.svg"
        ags.run_annotate(
            self.svg, self.csv, self.info, out, "SAMPLE01", meta)
        self.assertNotIn('data-bucket=', out.read_text())

    def test_bin_metadata_malformed_json_degrades(self):
        meta = _write(self.tmp, "bin_metadata.json", "{not json")
        out = self.tmp / "graph.svg"
        ags.run_annotate(
            self.svg, self.csv, self.info, out, "SAMPLE01", meta)
        self.assertNotIn('data-bucket=', out.read_text())

    def test_main_with_bin_metadata_flag(self):
        meta = _write(self.tmp, "bin_metadata.json", BIN_METADATA_JSON)
        out = self.tmp / "graph.svg"
        argv = [
            "annotate_graph_svg.py", "annotate",
            "--svg", str(self.svg), "--sentinels", str(self.csv),
            "--assembly-info", str(self.info), "--sample-id", "SAMPLE01",
            "--bin-metadata", str(meta), "--out", str(out),
        ]
        with patch.object(sys, "argv", argv):
            self.assertEqual(ags.main(), 0)
        self.assertIn('data-bucket="target"', out.read_text())


if __name__ == "__main__":
    unittest.main()
