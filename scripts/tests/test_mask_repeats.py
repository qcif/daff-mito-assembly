"""Unit tests for scripts/refdata/mask_repeats.py (task 68 §2)."""

import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
REFDATA = REPO_ROOT / "scripts" / "refdata"
sys.path.insert(0, str(REPO_ROOT / "bin"))
sys.path.insert(0, str(REFDATA))

_spec = importlib.util.spec_from_file_location(
    "mask_repeats", REFDATA / "mask_repeats.py")
mr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mr)

PARAMS = {'dust_level': 30, 'dust_window': 64, 'dust_linker': 1}


class MaskRepeatsTest(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = Path(self._tmp.name)

    def _run(self, records, dust_text, max_frac=0.6):
        fasta = self.dir / "in.fa"
        fasta.write_text("".join(
            f">{acc} desc\n{seq}\n" for acc, seq in records.items()))
        dust = self.dir / "in.dust"
        dust.write_text(dust_text)
        out = self.dir / "out.fa"
        with redirect_stderr(io.StringIO()):
            report = mr.mask_repeats(fasta, dust, out, max_frac, PARAMS)
        seqs = {}
        acc = None
        for line in out.read_text().splitlines():
            if line.startswith(">"):
                acc = line[1:].split()[0]
                seqs[acc] = ""
            else:
                seqs[acc] += line
        return seqs, report

    def test_interval_is_inclusive_and_converted(self):
        # dustmasker "2 - 4" covers bases 2, 3 and 4.
        self.assertEqual(
            self._parse(">a\n2 - 4\n"), {'a': [(2, 5)]})

    def _parse(self, text):
        path = self.dir / "p.dust"
        path.write_text(text)
        return mr.load_dust_intervals(path)

    def test_masks_interval_and_preserves_length(self):
        seqs, report = self._run(
            {'a': 'ACGTACGTAC'}, ">a\n2 - 4\n")
        self.assertEqual(seqs['a'], 'ACNNNCGTAC')
        self.assertEqual(report['bases_masked'], 3)

    def test_record_without_intervals_unchanged(self):
        seqs, report = self._run(
            {'a': 'ACGTACGTAC', 'b': 'TTTTGGGGCC'}, ">a\n0 - 1\n")
        self.assertEqual(seqs['b'], 'TTTTGGGGCC')
        self.assertEqual(report['genomes_with_hits'], 1)

    def test_blank_lines_ignored(self):
        self.assertEqual(
            self._parse("\n>a\n\n1 - 2\n"), {'a': [(1, 3)]})

    def test_guard_skips_mostly_masked_record(self):
        seqs, report = self._run(
            {'a': 'ACGTACGTAC', 'b': 'ACGTACGTAC'},
            ">a\n0 - 8\n>b\n0 - 1\n")
        self.assertNotIn('a', seqs)
        self.assertIn('b', seqs)
        self.assertEqual(report['genomes_skipped'], 1)

    def test_report_records_dust_parameters(self):
        _, report = self._run({'a': 'ACGT'}, "")
        self.assertEqual(report['masker'], 'dustmasker')
        self.assertEqual(report['dust_level'], 30)
        json.dumps(report)


if __name__ == "__main__":
    unittest.main()
