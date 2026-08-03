"""Unit tests for scripts/refdata/mask_panel.py — task 28 §4.

Build-time code, not a C1–C7 component, so CONSTITUTION rule 14's 100 %
branch-coverage requirement does not formally bind it. It is tested
anyway because a silent bug here produces a *plausible but wrong*
reference bundle, which then produces plausible but wrong results in
every downstream run — the hardest class of failure to notice
(task 28 §4).

The script is deliberately tool-free (minimap2 lives in the calling
build script), so these tests need no boundary mocking.

Cases:
  1. Interval merge: overlapping, adjacent, nested, disjoint, single,
     empty.
  2. Single interval masked mid-sequence; flanks untouched.
  3. Two overlapping intervals — union masked once, no double-count.
  4. Interval at sequence start / end — no off-by-one.
  5. Alignment shorter than --min-align-len — not masked.
  6. Genome above --max-masked-frac — skipped, warned, absent from
     output, recorded in the report.
  7. Genome with no plastid hits — passed through unchanged.
  8. Masked length == input length for every record (the invariant).
"""

import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest.mock import patch

REPO_ROOT = Path(__file__).resolve().parents[2]
BIN_DIR = REPO_ROOT / "bin"
SCRIPT = REPO_ROOT / "scripts" / "refdata" / "mask_panel.py"

# mask_panel.py imports the shared `intervals` module from bin/.
sys.path.insert(0, str(BIN_DIR))

_spec = importlib.util.spec_from_file_location("mask_panel", SCRIPT)
mp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mp)


def _write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def _fasta(records: dict) -> str:
    """Build a FASTA string from {accession: sequence}."""
    return "".join(
        f">{acc} synthetic record\n{mp.wrap(seq)}\n"
        for acc, seq in records.items()
    )


def _paf(rows: list) -> str:
    """
    Build a minimal PAF from (query, q_start, q_end) tuples.

    Only columns 1, 3 and 4 are read by load_mask_intervals; the rest
    are padded to a plausible 12-column PAF line.
    """
    return "".join(
        "\t".join([
            q, "0", str(qs), str(qe), "+",
            "pt_ref", "0", "0", str(qe - qs), str(qe - qs), "60",
        ]) + "\n"
        for q, qs, qe in rows
    )


class MaskPanelTestCase(unittest.TestCase):
    """Shared scratch dir + a runner that returns (report, records)."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def run_mask(self, records: dict, paf_rows: list, **kwargs):
        """Mask `records` against `paf_rows`; return (report, out, err)."""
        fasta = _write(self.tmp / "in.fa", _fasta(records))
        paf = _write(self.tmp / "hits.paf", _paf(paf_rows))
        out = self.tmp / "out.fa"
        params = {
            'min_align_len': mp.DEFAULT_MIN_ALIGN_LEN,
            'max_masked_frac': mp.DEFAULT_MAX_MASKED_FRAC,
        }
        params.update(kwargs)

        stderr = io.StringIO()
        with redirect_stderr(stderr):
            report = mp.mask_panel(
                fasta=fasta, paf=paf, out_fasta=out, **params)

        written = {
            mp.record_id(h): s
            for h, s in mp.iter_fasta_records(out)
        }
        return report, written, stderr.getvalue()


class TestMergeIntervals(MaskPanelTestCase):
    """Case 1 — the merge routine mask_panel shares with C2 and C3."""

    def test_merge_shapes(self):
        cases = [
            ("empty", [], [], 0),
            ("single", [(10, 20)], [(10, 20)], 10),
            ("disjoint", [(0, 10), (20, 30)], [(0, 10), (20, 30)], 20),
            ("overlapping", [(0, 15), (10, 25)], [(0, 25)], 25),
            ("adjacent", [(0, 10), (10, 20)], [(0, 20)], 20),
            ("nested", [(0, 100), (20, 30)], [(0, 100)], 100),
            ("unsorted", [(50, 60), (0, 10)], [(0, 10), (50, 60)], 20),
        ]
        for label, given, expect, length in cases:
            with self.subTest(label):
                self.assertEqual(mp.merge_intervals(given), expect)
                self.assertEqual(mp.merged_length(given), length)


class TestMasking(MaskPanelTestCase):
    """Cases 2–5, 7 — per-record masking behaviour."""

    def test_single_interval_mid_sequence(self):
        """Case 2: exactly the span becomes N; flanks byte-identical."""
        seq = "ACGT" * 250  # 1000 bp
        report, out, _ = self.run_mask(
            {"NC_000001.1": seq}, [("NC_000001.1", 300, 700)])

        masked = out["NC_000001.1"]
        self.assertEqual(masked[300:700], "N" * 400)
        self.assertEqual(masked[:300], seq[:300])
        self.assertEqual(masked[700:], seq[700:])
        self.assertEqual(report['bases_masked'], 400)

    def test_overlapping_intervals_counted_once(self):
        """Case 3: union masked once, no double-counting."""
        seq = "ACGT" * 250
        report, out, _ = self.run_mask(
            {"NC_000001.1": seq},
            [("NC_000001.1", 100, 400), ("NC_000001.1", 300, 600)],
        )

        self.assertEqual(out["NC_000001.1"].count("N"), 500)
        self.assertEqual(report['bases_masked'], 500)
        self.assertAlmostEqual(report['masked_fraction'], 0.5)

    def test_intervals_at_sequence_boundaries(self):
        """Case 4: start and end spans mask without off-by-one."""
        seq = "ACGT" * 250
        _, out, _ = self.run_mask(
            {"NC_000001.1": seq},
            [("NC_000001.1", 0, 250), ("NC_000001.1", 750, 1000)],
        )

        masked = out["NC_000001.1"]
        self.assertEqual(len(masked), 1000)
        self.assertEqual(masked[:250], "N" * 250)
        self.assertEqual(masked[250:750], seq[250:750])
        self.assertEqual(masked[750:], "N" * 250)

    def test_interval_clamped_to_sequence(self):
        """A PAF span running past the record end masks to the end."""
        seq = "ACGT" * 250
        report, out, _ = self.run_mask(
            {"NC_000001.1": seq}, [("NC_000001.1", 900, 1400)])

        self.assertEqual(len(out["NC_000001.1"]), 1000)
        self.assertEqual(out["NC_000001.1"][900:], "N" * 100)
        self.assertEqual(report['bases_masked'], 100)

    def test_short_alignment_not_masked(self):
        """Case 5: below --min-align-len is conserved gene, not NUPT."""
        seq = "ACGT" * 250
        report, out, _ = self.run_mask(
            {"NC_000001.1": seq},
            [("NC_000001.1", 300, 450)],  # 150 bp, floor is 200
            min_align_len=200,
        )

        self.assertEqual(out["NC_000001.1"], seq)
        self.assertEqual(report['bases_masked'], 0)
        self.assertEqual(report['genomes_with_plastid_hits'], 0)

    def test_alignment_exactly_at_floor_is_masked(self):
        """The floor is inclusive — 200 bp masks, 199 bp does not."""
        seq = "ACGT" * 250
        report, _, _ = self.run_mask(
            {"NC_000001.1": seq},
            [("NC_000001.1", 300, 500)],
            min_align_len=200,
        )
        self.assertEqual(report['bases_masked'], 200)

    def test_no_plastid_hits_passes_through(self):
        """Case 7: a genome with no hits is unchanged."""
        seq = "ACGTTGCA" * 125
        report, out, _ = self.run_mask({"NC_000002.1": seq}, [])

        self.assertEqual(out["NC_000002.1"], seq)
        self.assertEqual(report['bases_masked'], 0)
        self.assertEqual(report['genomes_written'], 1)
        self.assertEqual(report['genomes_skipped'], 0)


class TestMaxMaskedFractionGuard(MaskPanelTestCase):
    """Case 6 — the data-quality guard, not a tuning knob."""

    def test_genome_above_guard_is_skipped(self):
        seq = "ACGT" * 250
        report, out, err = self.run_mask(
            {"NC_000001.1": seq, "NC_000002.1": seq},
            [("NC_000001.1", 0, 900)],  # 90 %, above the 60 % guard
            max_masked_frac=0.60,
        )

        self.assertNotIn("NC_000001.1", out)
        self.assertIn("NC_000002.1", out)
        self.assertIn("skipping NC_000001.1", err)
        self.assertEqual(report['genomes_in'], 2)
        self.assertEqual(report['genomes_written'], 1)
        self.assertEqual(report['genomes_skipped'], 1)

        skipped = report['skipped_genomes']
        self.assertEqual([g['accession'] for g in skipped],
                         ["NC_000001.1"])
        self.assertAlmostEqual(skipped[0]['masked_fraction'], 0.9)

        # A skipped genome must not inflate the panel-wide totals.
        self.assertEqual(report['bases_written'], 1000)
        self.assertEqual(report['bases_masked'], 0)

    def test_genome_exactly_at_guard_is_kept(self):
        """The guard fires strictly above the fraction, not at it."""
        seq = "ACGT" * 250
        report, out, _ = self.run_mask(
            {"NC_000001.1": seq},
            [("NC_000001.1", 0, 600)],
            max_masked_frac=0.60,
        )

        self.assertIn("NC_000001.1", out)
        self.assertEqual(report['genomes_skipped'], 0)

    def test_top_masked_genomes_ranked(self):
        seq = "ACGT" * 250
        report, _, _ = self.run_mask(
            {"NC_000001.1": seq, "NC_000002.1": seq, "NC_000003.1": seq},
            [("NC_000001.1", 0, 200), ("NC_000002.1", 0, 500)],
        )

        ranked = [g['accession'] for g in report['top_masked_genomes']]
        self.assertEqual(
            ranked, ["NC_000002.1", "NC_000001.1", "NC_000003.1"])


class TestLengthInvariant(MaskPanelTestCase):
    """Case 8 — masking substitutes, never deletes."""

    def test_every_written_record_preserves_length(self):
        records = {
            "NC_000001.1": "ACGT" * 250,     # 1000
            "NC_000002.1": "TTGCA" * 300,    # 1500
            "NC_000003.1": "GGCCAT" * 100,   # 600
        }
        _, out, _ = self.run_mask(
            records,
            [
                ("NC_000001.1", 0, 400),
                ("NC_000002.1", 700, 1500),
                ("NC_000003.1", 100, 300),
                ("NC_000003.1", 250, 350),
            ],
        )

        self.assertEqual(set(out), set(records))
        for acc, seq in records.items():
            with self.subTest(acc):
                self.assertEqual(len(out[acc]), len(seq))

    def test_masked_bases_are_only_ever_added(self):
        """Non-masked positions keep their original base."""
        seq = "ACGTTGCAGGCC" * 100
        _, out, _ = self.run_mask(
            {"NC_000001.1": seq}, [("NC_000001.1", 200, 600)])

        masked = out["NC_000001.1"]
        for i, base in enumerate(masked):
            if not 200 <= i < 600:
                self.assertEqual(base, seq[i])


class TestCli(MaskPanelTestCase):
    """main() wiring — args parsed, report written, summary on stderr."""

    def test_main_writes_fasta_and_report(self):
        seq = "ACGT" * 250
        fasta = _write(self.tmp / "in.fa", _fasta({"NC_000001.1": seq}))
        paf = _write(self.tmp / "hits.paf", _paf([("NC_000001.1", 0, 400)]))
        out = self.tmp / "out.fa"
        report_path = self.tmp / "report.json"

        argv = [
            "mask_panel.py",
            "--fasta", str(fasta),
            "--paf", str(paf),
            "--out", str(out),
            "--report", str(report_path),
            "--min-align-len", "200",
            "--max-masked-frac", "0.6",
        ]
        stderr = io.StringIO()
        with patch.object(sys, "argv", argv):
            with redirect_stderr(stderr):
                mp.main()

        report = json.loads(report_path.read_text())
        self.assertEqual(report['min_align_len'], 200)
        self.assertEqual(report['max_masked_frac'], 0.6)
        self.assertEqual(report['bases_masked'], 400)
        self.assertIn("1 genomes written", stderr.getvalue())

        written = dict(mp.iter_fasta_records(out))
        self.assertEqual(len(next(iter(written.values()))), 1000)


if __name__ == "__main__":
    unittest.main()
