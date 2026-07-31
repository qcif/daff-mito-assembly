"""Unit tests for bin/plastid_canonicalise.py — spec §9 test-case matrix.

Fixtures are synthetic GFA strings built inline (no checked-in binary
GFAs). Sequence content is deterministic per edge name so the
identity assertions in cases #9/#10 are meaningful and reproducible.

Cases (spec/plastid-canonicalisation.md §9, numbered to match):
  1. 3-edge canonical.
  2. Upper-case `DP:f:` depth tag.
  3. Missing depth tags on all edges -> depth_tie.
  4. Missing depth tag on IR only -> lsc_ir_collision.
  5. LSC-IR collision (explicit depths).
  6. 1-edge resolved circle.
  7. 2-edge -> non_canonical.
  8. 4-edge -> non_canonical.
  9. Round-trip length invariant.
  10. path2 uses rc(SSC) at the correct offset.

Plus: parser edge cases (case-insensitive tag scan position, duplicate
edge names, `*` sequence, malformed depth float), tie-break branches
in edge assignment, zero-length-edge degenerate check, and the
construction invariant guard — for 100% branch coverage on the
parsing, classification and canonicalisation functions per
CONSTITUTION.md rule 14.
"""

import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.Seq import Seq

MODULE_PATH = (
    Path(__file__).resolve().parents[2] / "bin" / "plastid_canonicalise.py"
)

_spec = importlib.util.spec_from_file_location(
    "plastid_canonicalise", MODULE_PATH
)
pc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pc)


def _make_seq(name: str, length: int) -> str:
    """Deterministic, edge-distinct nucleotide sequence of given length."""
    bases = "ACGT"
    unit = "".join(bases[(ord(c) + i) % 4] for i, c in enumerate(name * 4))
    return (unit * (length // len(unit) + 1))[:length]


def _build_gfa(specs, tag="dp", type_char="f") -> tuple:
    """
    specs: list of (edge_name, length, depth_or_None).
    Returns (gfa_text, {edge_name: sequence}).
    """
    lines = ["H\tVN:Z:1.0"]
    sequences = {}
    for name, length, depth in specs:
        seq = _make_seq(name, length)
        sequences[name] = seq
        if depth is None:
            lines.append(f"S\t{name}\t{seq}")
        else:
            lines.append(f"S\t{name}\t{seq}\t{tag}:{type_char}:{depth}")
    lines.append("L\tedge_A\t+\tedge_B\t+\t0M")  # ignored link line
    return "\n".join(lines) + "\n", sequences


# Spec §9 case #1 fixture: LSC 90k depth 1.0, IR 25k depth 2.0, SSC 15k
# depth 1.0.
CASE1_SPECS = [
    ("edge_A", 90000, 1.0),
    ("edge_B", 25000, 2.0),
    ("edge_C", 15000, 1.0),
]


class TestCaseMatrix(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _write_gfa(self, text: str) -> Path:
        path = self.tmp / "assembly_graph.gfa"
        path.write_text(text)
        return path

    def test_case1_three_edge_canonical(self):
        gfa_text, _ = _build_gfa(CASE1_SPECS)
        gfa_path = self._write_gfa(gfa_text)
        outdir = self.tmp / "out"

        result = pc.canonicalise_plastid(gfa_path, outdir=outdir)

        self.assertEqual(result.branch, "canonical")
        self.assertEqual(result.lsc_edge, "edge_A")
        self.assertEqual(result.ir_edge, "edge_B")
        self.assertEqual(result.ssc_edge, "edge_C")
        self.assertTrue((outdir / "path1.fasta").exists())
        self.assertTrue((outdir / "path2.fasta").exists())

    def test_case2_uppercase_depth_tag(self):
        specs = [
            ("edge_A", 90000, 1.0),
            ("edge_B", 25000, None),
            ("edge_C", 15000, 1.0),
        ]
        gfa_text, _ = _build_gfa(specs)
        # Manually append the upper-case tag to edge_B's line.
        gfa_text = gfa_text.replace(
            f"S\tedge_B\t{_make_seq('edge_B', 25000)}",
            f"S\tedge_B\t{_make_seq('edge_B', 25000)}\tDP:f:2.0",
        )
        gfa_path = self._write_gfa(gfa_text)

        result = pc.canonicalise_plastid(gfa_path, outdir=self.tmp / "out")

        self.assertEqual(result.branch, "canonical")

    def test_case2b_integer_type_depth_tag(self):
        """Real Flye assembly_graph.gfa output uses dp:i:, not dp:f:."""
        specs = [
            ("edge_A", 90000, 1),
            ("edge_B", 25000, 2),
            ("edge_C", 15000, 1),
        ]
        gfa_text, _ = _build_gfa(specs, type_char="i")
        gfa_path = self._write_gfa(gfa_text)

        result = pc.canonicalise_plastid(gfa_path, outdir=self.tmp / "out")

        self.assertEqual(result.branch, "canonical")

    def test_case3_all_depths_missing_is_depth_tie(self):
        specs = [
            ("edge_A", 90000, None),
            ("edge_B", 25000, None),
            ("edge_C", 15000, None),
        ]
        gfa_text, _ = _build_gfa(specs)
        gfa_path = self._write_gfa(gfa_text)

        result = pc.canonicalise_plastid(gfa_path, outdir=self.tmp / "out")

        self.assertEqual(result.branch, "non_canonical")
        self.assertIn("depth_tie", result.non_canonical_reason)
        self.assertFalse((self.tmp / "out").exists())

    def test_case4_missing_depth_on_ir_only_collides(self):
        specs = [
            ("edge_A", 90000, 1.0),
            ("edge_B", 25000, None),
            ("edge_C", 15000, 1.0),
        ]
        gfa_text, _ = _build_gfa(specs)
        gfa_path = self._write_gfa(gfa_text)

        result = pc.canonicalise_plastid(gfa_path, outdir=self.tmp / "out")

        self.assertEqual(result.branch, "non_canonical")
        self.assertIn("lsc_ir_collision", result.non_canonical_reason)

    def test_case5_lsc_ir_collision_explicit(self):
        specs = [
            ("edge_A", 90000, 3.0),
            ("edge_B", 25000, 2.0),
            ("edge_C", 15000, 1.0),
        ]
        gfa_text, _ = _build_gfa(specs)
        gfa_path = self._write_gfa(gfa_text)

        result = pc.canonicalise_plastid(gfa_path, outdir=self.tmp / "out")

        self.assertEqual(result.branch, "non_canonical")
        self.assertEqual(
            result.non_canonical_reason,
            "lsc_ir_collision: longest and deepest are the same edge",
        )

    def test_case6_one_edge_resolved_circle(self):
        specs = [("edge_A", 160000, 50.0)]
        gfa_text, _ = _build_gfa(specs)
        gfa_path = self._write_gfa(gfa_text)
        outdir = self.tmp / "out"

        result = pc.canonicalise_plastid(gfa_path, outdir=outdir)

        self.assertEqual(result.branch, "resolved_circle")
        self.assertIsNone(result.lsc_edge)
        self.assertIsNone(result.non_canonical_reason)
        self.assertFalse(outdir.exists())

    def test_case7_two_edge_non_canonical(self):
        specs = [("edge_A", 90000, 1.0), ("edge_B", 25000, 2.0)]
        gfa_text, _ = _build_gfa(specs)
        gfa_path = self._write_gfa(gfa_text)

        result = pc.canonicalise_plastid(gfa_path, outdir=self.tmp / "out")

        self.assertEqual(result.branch, "non_canonical")
        self.assertEqual(result.non_canonical_reason, "edge_count:2")

    def test_case8_four_edge_non_canonical(self):
        specs = [
            ("edge_A", 90000, 1.0),
            ("edge_B", 25000, 2.0),
            ("edge_C", 15000, 1.0),
            ("edge_D", 5000, 1.0),
        ]
        gfa_text, _ = _build_gfa(specs)
        gfa_path = self._write_gfa(gfa_text)

        result = pc.canonicalise_plastid(gfa_path, outdir=self.tmp / "out")

        self.assertEqual(result.branch, "non_canonical")
        self.assertEqual(result.non_canonical_reason, "edge_count:4")

    def test_case9_round_trip_length_invariant(self):
        gfa_text, _ = _build_gfa(CASE1_SPECS)
        gfa_path = self._write_gfa(gfa_text)
        outdir = self.tmp / "out"

        result = pc.canonicalise_plastid(gfa_path, outdir=outdir)

        expected = 90000 + 2 * 25000 + 15000
        self.assertEqual(result.path1_len, expected)
        self.assertEqual(result.path2_len, expected)
        path1_seq = (outdir / "path1.fasta").read_text().splitlines()[1]
        path2_seq = (outdir / "path2.fasta").read_text().splitlines()[1]
        self.assertEqual(len(path1_seq), expected)
        self.assertEqual(len(path2_seq), expected)

    def test_case10_path2_uses_rc_ssc_at_correct_offset(self):
        gfa_text, sequences = _build_gfa(CASE1_SPECS)
        gfa_path = self._write_gfa(gfa_text)
        outdir = self.tmp / "out"

        pc.canonicalise_plastid(gfa_path, outdir=outdir)

        path2_seq = (outdir / "path2.fasta").read_text().splitlines()[1]
        ssc_rc = str(Seq(sequences["edge_C"]).reverse_complement())
        offset = 90000 + 25000
        self.assertEqual(path2_seq[offset:offset + 15000], ssc_rc)


class TestParsing(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _write_gfa(self, text: str) -> Path:
        path = self.tmp / "assembly_graph.gfa"
        path.write_text(text)
        return path

    def test_depth_tag_at_non_first_optional_position(self):
        seq = _make_seq("edge_A", 100)
        text = f"H\tVN:Z:1.0\nS\tedge_A\t{seq}\tRC:i:5\tdp:f:12.5\tLN:i:100\n"
        gfa_path = self._write_gfa(text)

        edges = pc.parse_gfa(gfa_path)

        self.assertEqual(len(edges), 1)
        self.assertEqual(edges[0].depth, 12.5)

    def test_duplicate_edge_name_raises(self):
        seq = _make_seq("edge_A", 100)
        text = f"S\tedge_A\t{seq}\tdp:f:1.0\nS\tedge_A\t{seq}\tdp:f:2.0\n"
        gfa_path = self._write_gfa(text)

        with self.assertRaises(ValueError):
            pc.parse_gfa(gfa_path)

    def test_star_sequence_raises(self):
        text = "S\tedge_A\t*\tdp:f:1.0\n"
        gfa_path = self._write_gfa(text)

        with self.assertRaises(ValueError):
            pc.parse_gfa(gfa_path)

    def test_malformed_depth_float_treated_as_zero(self):
        seq = _make_seq("edge_A", 100)
        text = f"S\tedge_A\t{seq}\tdp:f:1.2.3\n"
        gfa_path = self._write_gfa(text)

        edges = pc.parse_gfa(gfa_path)

        self.assertEqual(edges[0].depth, 0.0)

    def test_missing_gfa_file_raises_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            pc.canonicalise_plastid(self.tmp / "nonexistent.gfa")


class TestTieBreaks(unittest.TestCase):
    """Directly exercises the tie-break branches in _select_lsc/_select_ir."""

    def test_lsc_length_tie_broken_by_depth(self):
        edges = [
            pc.Edge("edge_A", "A" * 100, 1.0, 100),
            pc.Edge("edge_B", "C" * 100, 2.0, 100),
            pc.Edge("edge_C", "G" * 50, 1.0, 50),
        ]
        self.assertEqual(pc._select_lsc(edges).name, "edge_B")

    def test_lsc_full_tie_broken_by_name_ascending(self):
        edges = [
            pc.Edge("edge_Z", "A" * 100, 1.0, 100),
            pc.Edge("edge_A", "C" * 100, 1.0, 100),
        ]
        self.assertEqual(pc._select_lsc(edges).name, "edge_A")

    def test_ir_depth_tie_broken_by_length(self):
        edges = [
            pc.Edge("edge_A", "A" * 100, 2.0, 100),
            pc.Edge("edge_B", "C" * 50, 2.0, 50),
            pc.Edge("edge_C", "G" * 30, 1.0, 30),
        ]
        self.assertEqual(pc._select_ir(edges).name, "edge_A")

    def test_ir_full_tie_broken_by_name_descending(self):
        edges = [
            pc.Edge("edge_A", "A" * 100, 2.0, 100),
            pc.Edge("edge_Z", "C" * 100, 2.0, 100),
        ]
        self.assertEqual(pc._select_ir(edges).name, "edge_Z")


class TestDegenerateChecks(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _write_gfa(self, text: str) -> Path:
        path = self.tmp / "assembly_graph.gfa"
        path.write_text(text)
        return path

    def test_zero_length_edge_is_non_canonical(self):
        specs = [
            ("edge_A", 90000, 1.0),
            ("edge_B", 25000, 2.0),
            ("edge_C", 0, 1.0),
        ]
        # Zero-length sequence: bypass _make_seq (would divide fine at 0).
        gfa_text, _ = _build_gfa(specs)
        gfa_path = self._write_gfa(gfa_text)

        result = pc.canonicalise_plastid(gfa_path, outdir=self.tmp / "out")

        self.assertEqual(result.branch, "non_canonical")
        self.assertEqual(result.non_canonical_reason, "zero_length_edge")
        self.assertFalse((self.tmp / "out").exists())

    def test_invariant_failure_raises_value_error(self):
        """Construction invariant is defensive; force it via a broken RC."""
        gfa_text, _ = _build_gfa(CASE1_SPECS)
        gfa_path = self._write_gfa(gfa_text)

        class _BrokenSeq:
            def __init__(self, s):
                self.s = s

            def reverse_complement(self):
                return "TOO_SHORT"

        with patch.object(pc, "Seq", _BrokenSeq):
            with self.assertRaises(ValueError):
                pc.canonicalise_plastid(gfa_path, outdir=self.tmp / "out")


class TestCli(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_cli_writes_json_out(self):
        gfa_text, _ = _build_gfa(CASE1_SPECS)
        gfa_path = self.tmp / "assembly_graph.gfa"
        gfa_path.write_text(gfa_text)
        outdir = self.tmp / "out"
        json_out = self.tmp / "result.json"

        with patch(
            "sys.argv",
            ["plastid_canonicalise.py", str(gfa_path),
             "--outdir", str(outdir), "--json-out", str(json_out)],
        ):
            exit_code = pc.main()

        self.assertEqual(exit_code, 0)
        self.assertTrue(json_out.exists())

    def test_cli_non_canonical_summary(self):
        specs = [("edge_A", 90000, 1.0), ("edge_B", 25000, 2.0)]
        gfa_text, _ = _build_gfa(specs)
        gfa_path = self.tmp / "assembly_graph.gfa"
        gfa_path.write_text(gfa_text)

        with patch(
            "sys.argv",
            ["plastid_canonicalise.py", str(gfa_path)],
        ):
            exit_code = pc.main()

        self.assertEqual(exit_code, 0)

    def test_cli_resolved_circle_summary(self):
        specs = [("edge_A", 160000, 50.0)]
        gfa_text, _ = _build_gfa(specs)
        gfa_path = self.tmp / "assembly_graph.gfa"
        gfa_path.write_text(gfa_text)

        with patch(
            "sys.argv",
            ["plastid_canonicalise.py", str(gfa_path)],
        ):
            exit_code = pc.main()

        self.assertEqual(exit_code, 0)


if __name__ == "__main__":
    unittest.main()
