"""Unit tests for bin/coverage_gate.py — spec §2.1 decision branches.

Each test builds a synthetic gzipped FASTQ in a temp directory, invokes
coverage_gate.py as a subprocess, and asserts on the JSON outputs.
No Nextflow required.
"""

import gzip
import json
import subprocess
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "coverage_gate.py"
NOMINAL_SIZE = 17000  # animal_mt
MIN_COV = 30
MAX_COV = 300
READ_LEN = 100
SEED = 42


def make_fastq_gz(path: Path, n_reads: int) -> None:
    """Write a gzipped FASTQ with n_reads synthetic 100-bp reads."""
    with gzip.open(path, "wb") as fh:
        for i in range(n_reads):
            fh.write(f"@read_{i}\n".encode())
            fh.write((b"A" * READ_LEN) + b"\n")
            fh.write(b"+\n")
            fh.write((b"I" * READ_LEN) + b"\n")


def run_gate(tmp: Path, n_reads: int) -> tuple[dict, dict]:
    """Run coverage_gate.py with synthetic input; return (status, coverage)."""
    reads = tmp / "input.fastq.gz"
    make_fastq_gz(reads, n_reads)
    out_fastq = tmp / "gated.fastq.gz"
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--reads", str(reads),
            "--sample-id", "TEST",
            "--nominal-size", str(NOMINAL_SIZE),
            "--min-cov", str(MIN_COV),
            "--max-cov", str(MAX_COV),
            "--seed", str(SEED),
            "--out-fastq", str(out_fastq),
        ],
        capture_output=True,
        text=True,
        cwd=tmp,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"coverage_gate.py exited {result.returncode}:\n{result.stderr}"
        )
    status = json.loads((tmp / "sample_status.json").read_text())
    coverage = json.loads((tmp / "coverage.json").read_text())
    return status, coverage


class TestCoverageGate(unittest.TestCase):

    def setUp(self):
        import tempfile
        self._tmpdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmpdir.name)

    def tearDown(self):
        self._tmpdir.cleanup()

    def test_passthrough(self):
        """Case 1: ~150× input → passthrough, status=ok."""
        # 150 × 17000 / 100 = 25500 reads → ~150× on 17 kb nominal
        n_reads = 25500
        status, coverage = run_gate(self.tmp, n_reads)
        self.assertEqual(status["status"], "ok")
        self.assertFalse(coverage["subsampled"])
        self.assertEqual(coverage["fraction"], 1.0)
        # estimated_cov should be between MIN and MAX
        self.assertGreaterEqual(status["estimated_cov"], MIN_COV)
        self.assertLessEqual(status["estimated_cov"], MAX_COV)

    def test_low_coverage(self):
        """Case 2: 100 kb input (≈6×) → low_coverage."""
        # 100 reads × 100 bp = 10 000 bases → ≈0.59× on 17 kb nominal
        n_reads = 100
        status, coverage = run_gate(self.tmp, n_reads)
        self.assertEqual(status["status"], "low_coverage")
        self.assertFalse(coverage["subsampled"])
        self.assertLess(status["estimated_cov"], MIN_COV)

    def test_subsample(self):
        """Case 3: 10 Mb input (≈588×) → subsampled to ~300×."""
        # 100 000 reads × 100 bp = 10 000 000 bases → ≈588× on 17 kb nominal
        n_reads = 100000
        status, coverage = run_gate(self.tmp, n_reads)
        self.assertEqual(status["status"], "ok")
        self.assertTrue(coverage["subsampled"])
        self.assertLess(coverage["fraction"], 1.0)
        # post_subsample_cov should be within ±10% of MAX_COV
        self.assertAlmostEqual(
            coverage["post_subsample_cov"], MAX_COV, delta=MAX_COV * 0.1
        )

    def test_empty_input(self):
        """Case 4: 0 reads → low_coverage."""
        n_reads = 0
        status, coverage = run_gate(self.tmp, n_reads)
        self.assertEqual(status["status"], "low_coverage")
        self.assertEqual(status["total_recruited_bases"], 0)
        self.assertEqual(status["estimated_cov"], 0.0)


if __name__ == "__main__":
    unittest.main()
