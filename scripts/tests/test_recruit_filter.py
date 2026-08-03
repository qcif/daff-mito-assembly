"""Unit tests for bin/recruit_filter.py — spec §2 stage 5, §9 item 1.

C-component runtime code, so CONSTITUTION rule 14's 100 % branch
coverage of the selection path applies. Synthetic inline PAF only; the
script is deliberately tool-free (minimap2 runs in the module), so no
boundary mocking is needed.

Cases:
  1. Merged extent, not summed extent — overlapping blocks counted once.
  2. Multi-block read: disjoint blocks summed.
  3. min_aligned_frac admits/rejects at the boundary.
  4. min_aligned_bp admits/rejects at the boundary.
  5. Both floors apply — failing either rejects.
  6. Defaults of 0/0 recruit every aligned read.
  7. A read absent from the PAF is never recruited.
  8. Output is sorted and newline-terminated; stats JSON matches.
  9. Zero-length read record does not divide by zero.
"""

import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

REPO_ROOT = Path(__file__).resolve().parents[2]
BIN_DIR = REPO_ROOT / "bin"
SCRIPT = BIN_DIR / "recruit_filter.py"

# recruit_filter.py imports the shared `intervals` module from bin/.
sys.path.insert(0, str(BIN_DIR))

_spec = importlib.util.spec_from_file_location("recruit_filter", SCRIPT)
rf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rf)


def _paf(rows) -> str:
    """Build a PAF from (query, qlen, qstart, qend) tuples."""
    return "".join(
        "\t".join([
            q, str(qlen), str(qs), str(qe), "+",
            "panel_ref", "500000", "0", str(qe - qs),
            str(qe - qs), "60",
        ]) + "\n"
        for q, qlen, qs, qe in rows
    )


class RecruitFilterTestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def load(self, rows):
        paf = self.tmp / "hits.paf"
        paf.write_text(_paf(rows))
        return rf.read_alignments(paf)


class TestMergedExtent(RecruitFilterTestCase):
    def test_overlapping_blocks_counted_once(self):
        """Case 1: the metric is merged, not summed."""
        aln = self.load([
            ("read_a", 10000, 0, 3000),
            ("read_a", 10000, 2000, 5000),
        ])
        self.assertEqual(aln["read_a"], (10000, 5000))

    def test_disjoint_blocks_summed(self):
        """Case 2: separate blocks both count."""
        aln = self.load([
            ("read_a", 10000, 0, 1000),
            ("read_a", 10000, 8000, 9000),
        ])
        self.assertEqual(aln["read_a"], (10000, 2000))

    def test_blank_lines_ignored(self):
        paf = self.tmp / "hits.paf"
        paf.write_text(_paf([("read_a", 100, 0, 50)]) + "\n\n")
        self.assertEqual(rf.read_alignments(paf), {"read_a": (100, 50)})


class TestThresholds(RecruitFilterTestCase):
    def test_min_aligned_frac_boundary(self):
        """Case 3: the floor is inclusive."""
        aln = {"exact": (1000, 300), "below": (1000, 299)}
        self.assertEqual(rf.select(aln, 0.3, 0), ["exact"])

    def test_min_aligned_bp_boundary(self):
        """Case 4: the floor is inclusive."""
        aln = {"exact": (10000, 1000), "below": (10000, 999)}
        self.assertEqual(rf.select(aln, 0.0, 1000), ["exact"])

    def test_both_floors_apply(self):
        """Case 5: clearing one but not the other is a reject."""
        aln = {
            "both":     (10000, 5000),   # frac 0.50, bp 5000 → keep
            "frac_only": (2000, 1800),   # frac 0.90, bp 1800 → bp fails
            "bp_only":  (100000, 4000),  # frac 0.04, bp 4000 → frac fails
        }
        self.assertEqual(rf.select(aln, 0.3, 2000), ["both"])

    def test_defaults_recruit_every_aligned_read(self):
        """Case 6: 0/0 is a pass-through — the shipped default."""
        aln = {"a": (10000, 1), "b": (10000, 9999)}
        self.assertEqual(rf.select(aln, 0.0, 0), ["a", "b"])

    def test_unaligned_read_never_recruited(self):
        """Case 7: absence from the PAF means absence from the output."""
        aln = self.load([("read_a", 10000, 0, 5000)])
        self.assertEqual(rf.select(aln, 0.0, 0), ["read_a"])
        self.assertNotIn("read_b", aln)

    def test_zero_length_record_does_not_divide_by_zero(self):
        """Case 9: a malformed 0-length record must not crash."""
        aln = {"zero": (0, 0)}
        self.assertEqual(rf.select(aln, 0.5, 0), ["zero"])
        self.assertEqual(rf.select(aln, 0.5, 1), [])

    def test_output_sorted(self):
        aln = {"c": (100, 100), "a": (100, 100), "b": (100, 100)}
        self.assertEqual(rf.select(aln, 0.0, 0), ["a", "b", "c"])


class TestCli(RecruitFilterTestCase):
    def _run(self, rows, extra=()):
        paf = self.tmp / "hits.paf"
        paf.write_text(_paf(rows))
        out = self.tmp / "ids.txt"
        stats = self.tmp / "stats.json"
        argv = [
            "recruit_filter.py",
            "--paf", str(paf),
            "--out", str(out),
            "--stats", str(stats),
            *extra,
        ]
        stderr = io.StringIO()
        with patch.object(sys, "argv", argv), redirect_stderr(stderr):
            rf.main()
        return out, stats, stderr.getvalue()

    def test_writes_ids_and_stats(self):
        """Case 8: file output + stats sidecar agree."""
        out, stats, err = self._run(
            [("read_a", 10000, 0, 5000), ("read_b", 10000, 0, 100)],
            extra=["--min-aligned-frac", "0.3"],
        )
        self.assertEqual(out.read_text(), "read_a\n")
        self.assertIn("1 of 2 aligned reads recruited", err)

        recorded = json.loads(stats.read_text())
        self.assertEqual(recorded['reads_aligned'], 2)
        self.assertEqual(recorded['reads_recruited'], 1)
        self.assertEqual(recorded['min_aligned_frac'], 0.3)
        self.assertEqual(recorded['min_aligned_bp'], 0)

    def test_defaults_pass_everything_through(self):
        out, stats, _ = self._run(
            [("read_a", 10000, 0, 1), ("read_b", 10000, 0, 9999)])
        self.assertEqual(out.read_text(), "read_a\nread_b\n")
        self.assertEqual(json.loads(stats.read_text())['reads_recruited'], 2)

    def test_stdout_output_and_no_stats(self):
        """--out '-' writes to stdout; --stats omitted writes no file."""
        paf = self.tmp / "hits.paf"
        paf.write_text(_paf([("read_a", 100, 0, 100)]))
        argv = ["recruit_filter.py", "--paf", str(paf)]
        buf, stderr = io.StringIO(), io.StringIO()
        with patch.object(sys, "argv", argv), \
                redirect_stdout(buf), redirect_stderr(stderr):
            rf.main()
        self.assertEqual(buf.getvalue(), "read_a\n")
        self.assertFalse((self.tmp / "stats.json").exists())

    def test_empty_paf_yields_empty_selection(self):
        paf = self.tmp / "hits.paf"
        paf.write_text("")
        out = self.tmp / "ids.txt"
        argv = ["recruit_filter.py", "--paf", str(paf), "--out", str(out)]
        stderr = io.StringIO()
        with patch.object(sys, "argv", argv), redirect_stderr(stderr):
            rf.main()
        self.assertEqual(out.read_text(), "")
        self.assertIn("0 of 0", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
