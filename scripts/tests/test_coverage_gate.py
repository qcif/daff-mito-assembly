"""Unit tests for bin/coverage_gate.py — spec §2.1 decision branches.

Both `TestCoverageGate` and `TestSiblingSplit` load coverage_gate.py
**in-process** (`importlib`, module registered in `sys.modules` so
`patch()` can address it by name) rather than invoking it as a
subprocess. This is what makes the module's branches show up in the
`coverage` report at all: `scripts/pytest.sh` traces only the parent
process, so a subprocess-invoked script contributes zero measured
coverage no matter how many assertions its tests make (task 27).

Two notes on the mocking, easy to miss:

* `coverage_gate.py` does `import subprocess`, **not**
  `from subprocess import run`, so `patch('coverage_gate.subprocess.run')`
  patches `run` on the shared module object for the duration of the
  `with` block. Harmless here (single-threaded, nothing else in the
  module shells out) but non-obvious.
* `coverage_gate.py` writes `sample_status.json` and `coverage.json` to
  **CWD**. `contextlib.chdir` (3.11+; the image is 3.12) replaces what
  `cwd=tmp` gave the old subprocess form.

Residual risk — faking the seqkit/seqtk/minimap2 boundary means these
tests do **not** catch:

1. A renamed `sum_len` column in real `seqkit stats -T` output (the
   reordered-columns case below pins column *position* independence,
   not the name itself).
2. Real-tool edge behaviour on real input — most importantly whether a
   future seqkit version emits a data row for a zero-byte gzip.
3. The unquoted f-string interpolation of `args.reads` / `args.out_fastq`
   into the `bash -c` line, which breaks on any path containing a space
   or shell metacharacter. The bash-command-shape case below pins the
   string coverage_gate.py *builds*, not that a shell accepts it.

Per CONSTITUTION.md rule 15 these belong to the nightly
`-profile integration` surface, which runs the real
`neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5`
container — see tests/integration/assertions.sh.
"""

import contextlib
import gzip
import importlib.util
import json
import random
import re
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
SCRIPT = BIN_DIR / "coverage_gate.py"

# coverage_gate.py imports the shared `intervals` module as a bin/
# sibling, relying on Nextflow's bin/ staging at runtime.
sys.path.insert(0, str(BIN_DIR))

NOMINAL_SIZE = 17000  # animal_mt
MIN_COV = 30
MAX_COV = 300
READ_LEN = 100
SEED = 42


def _load_module():
    """Import coverage_gate.py in-process so patch() can address it."""
    spec = importlib.util.spec_from_file_location("coverage_gate", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["coverage_gate"] = module
    spec.loader.exec_module(module)
    return module


cg = _load_module()


def make_fastq_gz(path: Path, n_reads: int) -> None:
    """Write a gzipped FASTQ with n_reads synthetic 100-bp reads."""
    with gzip.open(path, "wb") as fh:
        for i in range(n_reads):
            fh.write(f"@read_{i}\n".encode())
            fh.write((b"A" * READ_LEN) + b"\n")
            fh.write(b"+\n")
            fh.write((b"I" * READ_LEN) + b"\n")


def _stats_table(path: Path) -> str:
    """Render a real `seqkit stats -T` table from the real gzip at path.

    Data-driven, not canned: the row is computed from what's actually
    in the file, so the coverage arithmetic the tests assert on is a
    genuine computation rather than a tautology. Emits a data row with
    sum_len=0 for an empty gzip, matching real seqkit's behaviour —
    total_bases() indexes out[1] unconditionally, so a header-only
    fake would IndexError rather than pass test_empty_input for the
    right reason.
    """
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    n_seqs = len(lines) // 4
    sum_len = sum(len(lines[i * 4 + 1].rstrip("\n"))
                  for i in range(n_seqs))
    header = "file\tformat\ttype\tnum_seqs\tsum_len"
    row = f"{path}\tFASTQ\tDNA\t{n_seqs}\t{sum_len}"
    return header + "\n" + row + "\n"


_SEQTK_CMD_RE = re.compile(
    r"seqtk sample -s (\d+) (\S+) ([0-9.eE+-]+) \| gzip > (\S+)")


def _fake_seqtk(shell_cmd: str) -> None:
    """Emulate `seqtk sample -s S IN F | gzip > OUT`.

    Bernoulli-selects records with random.Random(seed) — real seqtk's
    sampling semantics, not an exact count — and writes a real gzip so
    a subsequent _stats_table() call re-reads genuine content.
    """
    seed, reads_path, frac, out_path = _SEQTK_CMD_RE.match(shell_cmd).groups()
    rnd = random.Random(int(seed))
    frac = float(frac)
    with gzip.open(reads_path, "rt") as fh:
        lines = fh.readlines()
    records = [lines[i:i + 4] for i in range(0, len(lines), 4)]
    kept = [record for record in records if rnd.random() < frac]
    with gzip.open(out_path, "wt") as fh:
        for record in kept:
            fh.writelines(record)


def _fake_run(cmd, **kwargs):
    """Default fake for coverage_gate.subprocess.run.

    Raises on anything unrecognised so a future third tool call fails
    loudly rather than silently returning a MagicMock that satisfies
    every assertion.
    """
    if cmd[0] == "seqkit":
        return subprocess.CompletedProcess(
            cmd, 0, stdout=_stats_table(Path(cmd[3])))
    if cmd[0] == "bash":
        _fake_seqtk(cmd[2])
        return subprocess.CompletedProcess(cmd, 0)
    raise AssertionError(f"unexpected command: {cmd}")


def _run_gate(tmp: Path, n_reads: int, fake_run=None,
              **argv_overrides) -> tuple:
    """Drive coverage_gate.py's main() in-process against synthetic input.

    Returns (rc, status, coverage). `argv_overrides` lets individual
    tests override any CLI argument by its underscored name (e.g.
    `nominal_size="0"`).
    """
    reads = tmp / "input.fastq.gz"
    make_fastq_gz(reads, n_reads)
    out_fastq = tmp / "gated.fastq.gz"
    argv_map = {
        "reads": str(reads),
        "sample_id": "TEST",
        # animal_mt has no sibling panel, so the split is skipped and
        # no reference bundle is read — see SIBLING_PANELS.
        "assembly_target": "animal_mt",
        "ref_dir": str(tmp),
        "nominal_size": str(NOMINAL_SIZE),
        "min_cov": str(MIN_COV),
        "max_cov": str(MAX_COV),
        "seed": str(SEED),
        "out_fastq": str(out_fastq),
    }
    argv_map.update(argv_overrides)
    argv = ["coverage_gate.py"]
    for key, value in argv_map.items():
        argv += [f"--{key.replace('_', '-')}", value]

    with contextlib.chdir(tmp), \
            patch.object(sys, "argv", argv), \
            patch("coverage_gate.subprocess.run",
                  side_effect=fake_run or _fake_run):
        rc = cg.main()

    status = json.loads((tmp / "sample_status.json").read_text())
    coverage = json.loads((tmp / "coverage.json").read_text())
    return rc, status, coverage


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
        rc, status, coverage = _run_gate(self.tmp, n_reads)
        self.assertEqual(rc, 0)
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
        rc, status, coverage = _run_gate(self.tmp, n_reads)
        self.assertEqual(rc, 0)
        self.assertEqual(status["status"], "low_coverage")
        self.assertFalse(coverage["subsampled"])
        self.assertLess(status["estimated_cov"], MIN_COV)

    def test_subsample(self):
        """Case 3: 10 Mb input (≈588×) → subsampled to ~300×."""
        # 100 000 reads × 100 bp = 10 000 000 bases → ≈588× on 17 kb nominal
        n_reads = 100000
        rc, status, coverage = _run_gate(self.tmp, n_reads)
        self.assertEqual(rc, 0)
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
        rc, status, coverage = _run_gate(self.tmp, n_reads)
        self.assertEqual(rc, 0)
        self.assertEqual(status["status"], "low_coverage")
        self.assertEqual(status["total_recruited_bases"], 0)
        self.assertEqual(status["estimated_cov"], 0.0)

    def test_nominal_size_zero(self):
        """Case 5: --nominal-size 0 takes the false arm of the ternary."""
        rc, status, coverage = _run_gate(
            self.tmp, 100, nominal_size="0")
        self.assertEqual(rc, 0)
        self.assertEqual(status["estimated_cov"], 0.0)
        self.assertEqual(coverage["pre_subsample_cov"], 0.0)

    def test_exactly_max_cov_is_passthrough(self):
        """Case 6: est == max_cov exactly — pins `>` vs `>=` at the top."""
        # 300 × 17000 / 100 = 51 000 reads → exactly 300×
        rc, status, coverage = _run_gate(self.tmp, 51000)
        self.assertEqual(rc, 0)
        self.assertEqual(status["estimated_cov"], MAX_COV)
        self.assertEqual(status["status"], "ok")
        self.assertFalse(coverage["subsampled"])

    def test_exactly_min_cov_is_passthrough(self):
        """Case 7: est == min_cov exactly — pins `<` vs `<=` at the floor."""
        # 30 × 17000 / 100 = 5100 reads → exactly 30×
        rc, status, coverage = _run_gate(self.tmp, 5100)
        self.assertEqual(rc, 0)
        self.assertEqual(status["estimated_cov"], MIN_COV)
        self.assertEqual(status["status"], "ok")

    def test_stats_table_column_reordered(self):
        """Case 8: total_bases() finds sum_len by name, not position 5.

        Newer seqkit adds Q1/Q2 trailing columns; a hard-coded column-5
        read would silently pick up the wrong value instead of erroring.
        """
        def fake_run(cmd, **kwargs):
            self.assertEqual(cmd[0], "seqkit")
            return subprocess.CompletedProcess(
                cmd, 0,
                stdout="file\tformat\ttype\tsum_len\tnum_seqs\tQ1\tQ2\n"
                       "in\tFASTQ\tDNA\t12345\t2\t10\t20\n")
        with patch("coverage_gate.subprocess.run", side_effect=fake_run):
            self.assertEqual(
                cg.total_bases(self.tmp / "in.fastq.gz"), 12345)

    def test_subsample_bash_command_shape(self):
        """Case 9: the `bash -c` string mocking would otherwise erase.

        Asserts `seqtk sample -s 42`, the reads path, the fraction, and
        `| gzip > <out>` — the CLI contract, not just its effect.
        """
        recorded = {}

        def fake_run(cmd, **kwargs):
            if cmd[0] == "seqkit":
                return subprocess.CompletedProcess(
                    cmd, 0, stdout=_stats_table(Path(cmd[3])))
            if cmd[0] == "bash":
                recorded["shell_cmd"] = cmd[2]
                _fake_seqtk(cmd[2])
                return subprocess.CompletedProcess(cmd, 0)
            raise AssertionError(f"unexpected command: {cmd}")

        rc, status, coverage = _run_gate(
            self.tmp, 100000, fake_run=fake_run)
        self.assertEqual(rc, 0)
        self.assertTrue(coverage["subsampled"])
        shell_cmd = recorded["shell_cmd"]
        self.assertIn(f"seqtk sample -s {SEED}", shell_cmd)
        self.assertIn(str(self.tmp / "input.fastq.gz"), shell_cmd)
        self.assertIn("| gzip >", shell_cmd)
        self.assertIn(str(self.tmp / "gated.fastq.gz"), shell_cmd)

    def test_tool_failure_propagates_nonzero_exit(self):
        """Case 10: a real tool crash is not swallowed into "ok".

        Non-zero exits are reserved for unexpected tool failures — the
        module docstring's contract Nextflow's errorStrategy='ignore'
        depends on — so a CalledProcessError must propagate, not be
        caught and turned into a coverage decision.
        """
        def boom(cmd, **kwargs):
            raise subprocess.CalledProcessError(1, cmd)

        with self.assertRaises(subprocess.CalledProcessError):
            _run_gate(self.tmp, 100, fake_run=boom)


def paf_line(read_id, read_len, start, end, target="ref") -> str:
    """One minimal PAF record — only fields 1-4 are read by C2."""
    return (
        f"{read_id}\t{read_len}\t{start}\t{end}\t+\t{target}\t"
        f"1000\t0\t{end - start}\t{end - start}\t{end - start}\t60\n"
    )


class TestSiblingSplit(unittest.TestCase):
    """Task 25 §3.2 — sibling-organelle split of the recruited pool."""

    def setUp(self):
        import tempfile
        self._tmpdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmpdir.name)

    def tearDown(self):
        self._tmpdir.cleanup()

    # -- merge_intervals -------------------------------------------------

    def test_merge_intervals_shapes(self):
        """Overlapping, adjacent, nested, disjoint, single, empty."""
        cases = [
            ([], []),
            ([(0, 10)], [(0, 10)]),
            ([(0, 10), (5, 20)], [(0, 20)]),          # overlapping
            ([(0, 10), (10, 20)], [(0, 20)]),         # adjacent
            ([(0, 100), (10, 20)], [(0, 100)]),       # nested
            ([(0, 10), (20, 30)], [(0, 10), (20, 30)]),  # disjoint
        ]
        for intervals, expected in cases:
            with self.subTest(intervals=intervals):
                self.assertEqual(cg.merge_intervals(intervals), expected)

    def test_merge_intervals_unsorted(self):
        """Input order does not matter."""
        self.assertEqual(
            cg.merge_intervals([(20, 30), (0, 10), (5, 25)]), [(0, 30)])

    # -- split_by_panel --------------------------------------------------

    def test_split_target_wins(self):
        """A read aligning better to the declared panel is target-assigned."""
        hits = {
            "plant_mt": {"r1": (1000, 800)},
            "plant_pt": {"r1": (1000, 100)},
        }
        split = cg.split_by_panel(hits, "plant_mt")
        self.assertEqual(split["target_assigned_bases"], 1000)
        self.assertEqual(split["sibling_assigned_bases"], 0)
        self.assertEqual(split["target_merged_aligned_bases"], 800)
        self.assertEqual(split["reads_target_assigned"], 1)

    def test_split_sibling_wins(self):
        """Plastid carry-over is charged to the sibling, not the target."""
        hits = {
            "plant_mt": {"r1": (1000, 50)},
            "plant_pt": {"r1": (1000, 1000)},
        }
        split = cg.split_by_panel(hits, "plant_mt")
        self.assertEqual(split["target_assigned_bases"], 0)
        self.assertEqual(split["sibling_assigned_bases"], 1000)
        self.assertEqual(split["reads_sibling_assigned"], 1)

    def test_split_tie_falls_to_sibling(self):
        """A read equally explained by both is not target depth (rule 18)."""
        hits = {
            "plant_mt": {"r1": (1000, 500)},
            "plant_pt": {"r1": (1000, 500)},
        }
        split = cg.split_by_panel(hits, "plant_mt")
        self.assertEqual(split["target_assigned_bases"], 0)
        self.assertEqual(split["sibling_assigned_bases"], 1000)

    def test_split_read_absent_from_declared_panel(self):
        """Read length is recovered from the sibling hit when unmapped."""
        hits = {
            "plant_mt": {},
            "plant_pt": {"r1": (1234, 900)},
        }
        split = cg.split_by_panel(hits, "plant_mt")
        self.assertEqual(split["sibling_assigned_bases"], 1234)
        self.assertEqual(split["target_assigned_bases"], 0)

    # -- map_to_panel ----------------------------------------------------

    def test_map_to_panel_merges_blocks(self):
        """Multiple alignment blocks for one read collapse to their union."""
        paf = self.tmp / "out.paf"

        def fake_run(cmd, stdout=None, **kwargs):
            stdout.write(paf_line("r1", 1000, 0, 300))
            stdout.write(paf_line("r1", 1000, 200, 500))   # overlaps
            stdout.write(paf_line("r2", 500, 0, 100))
            return subprocess.CompletedProcess(cmd, 0)

        with patch("coverage_gate.subprocess.run", side_effect=fake_run):
            hits = cg.map_to_panel(
                Path("reads.fq.gz"), Path("panel.mmi"), paf, 1)

        self.assertEqual(hits["r1"], (1000, 500))   # 0-500 merged
        self.assertEqual(hits["r2"], (500, 100))

    # -- estimate_target_bases fallbacks ---------------------------------

    def test_no_sibling_panel_defined(self):
        """animal_mt has no sibling — split skipped, no mapping attempted."""
        target_bases, scored, split = cg.estimate_target_bases(
            Path("reads.fq.gz"), self.tmp, "animal_mt", 1)
        self.assertIsNone(target_bases)
        self.assertEqual(scored, [])
        self.assertEqual(split, {})

    def test_missing_sibling_index_falls_back(self):
        """Sibling .mmi absent → warn, no split, sample continues."""
        (self.tmp / "plant_mt.mmi").write_bytes(b"")   # declared present
        target_bases, scored, split = cg.estimate_target_bases(
            Path("reads.fq.gz"), self.tmp, "plant_mt", 1)
        self.assertIsNone(target_bases)
        self.assertEqual(scored, [])

    def test_missing_declared_index_falls_back(self):
        """Declared .mmi absent → same fallback, never a hard failure."""
        (self.tmp / "plant_pt.mmi").write_bytes(b"")
        target_bases, scored, _ = cg.estimate_target_bases(
            Path("reads.fq.gz"), self.tmp, "plant_mt", 1)
        self.assertIsNone(target_bases)
        self.assertEqual(scored, [])

    def test_minimap2_failure_falls_back(self):
        """minimap2 crash → warn, whole-pool estimate, no exception."""
        for panel in ("plant_mt", "plant_pt"):
            (self.tmp / f"{panel}.mmi").write_bytes(b"")

        def boom(cmd, **kwargs):
            raise subprocess.CalledProcessError(1, cmd)

        with patch("coverage_gate.subprocess.run", side_effect=boom):
            target_bases, scored, _ = cg.estimate_target_bases(
                Path("reads.fq.gz"), self.tmp, "plant_mt", 1)
        self.assertIsNone(target_bases)
        self.assertEqual(scored, [])

    def test_split_applied_when_both_panels_present(self):
        """Both panels load → estimate is target-assigned bases only."""
        for panel in ("plant_mt", "plant_pt"):
            (self.tmp / f"{panel}.mmi").write_bytes(b"")

        def fake_run(cmd, stdout=None, **kwargs):
            # cmd is the minimap2 argv; the .mmi is the second-last arg.
            panel = Path(cmd[-2]).stem
            if panel == "plant_mt":
                stdout.write(paf_line("mito", 1000, 0, 900))
                stdout.write(paf_line("plastid", 2000, 0, 40))
            else:
                stdout.write(paf_line("plastid", 2000, 0, 2000))
            return subprocess.CompletedProcess(cmd, 0)

        with patch("coverage_gate.subprocess.run", side_effect=fake_run):
            target_bases, scored, split = cg.estimate_target_bases(
                Path("reads.fq.gz"), self.tmp, "plant_mt", 4)

        self.assertEqual(target_bases, 1000)          # the mito read only
        self.assertEqual(scored, ["plant_pt"])
        self.assertEqual(split["sibling_assigned_bases"], 2000)

    # -- main() end to end, tools faked ----------------------------------

    def _run_main(self, tmp, fake_run, target="plant_mt", nominal=400000):
        reads = tmp / "input.fastq.gz"
        make_fastq_gz(reads, 10)
        argv = [
            "coverage_gate.py",
            "--reads", str(reads),
            "--sample-id", "TEST",
            "--assembly-target", target,
            "--ref-dir", str(tmp),
            "--nominal-size", str(nominal),
            "--min-cov", str(MIN_COV),
            "--max-cov", str(MAX_COV),
            "--seed", str(SEED),
            "--out-fastq", str(tmp / "gated.fastq.gz"),
        ]
        import contextlib
        with contextlib.chdir(tmp), \
                patch.object(sys, "argv", argv), \
                patch("coverage_gate.subprocess.run", side_effect=fake_run):
            rc = cg.main()
        return (
            rc,
            json.loads((tmp / "sample_status.json").read_text()),
            json.loads((tmp / "coverage.json").read_text()),
        )

    def test_main_gates_on_target_bases(self):
        """The INT-PLANT-01-mt shape: plastid-dominated pool soft-fails.

        20 Mb recruited would read as 50× on a 400 kb nominal mitogenome
        and pass. Only 4 Mb is mitochondrial → 10× → low_coverage.
        """
        for panel in ("plant_mt", "plant_pt"):
            (self.tmp / f"{panel}.mmi").write_bytes(b"")

        def fake_run(cmd, stdout=None, **kwargs):
            if cmd[0] == "seqkit":
                return subprocess.CompletedProcess(
                    cmd, 0,
                    stdout="file\tformat\ttype\tnum_seqs\tsum_len\n"
                           "in\tFASTQ\tDNA\t2\t20000000\n")
            panel = Path(cmd[-2]).stem
            if panel == "plant_mt":
                stdout.write(paf_line("mito", 4000000, 0, 3900000))
                stdout.write(paf_line("plastid", 16000000, 0, 100000))
            else:
                stdout.write(paf_line("plastid", 16000000, 0, 15000000))
            return subprocess.CompletedProcess(cmd, 0)

        rc, status, coverage = self._run_main(self.tmp, fake_run)

        self.assertEqual(rc, 0)
        self.assertEqual(status["status"], "low_coverage")
        self.assertEqual(status["coverage_basis"], cg.BASIS_TARGET_ASSIGNED)
        self.assertEqual(status["estimated_cov"], 10.0)
        self.assertEqual(status["total_recruited_bases"], 20000000)
        self.assertEqual(status["target_assigned_bases"], 4000000)
        self.assertEqual(status["sibling_organelle_fraction"], 0.8)
        self.assertEqual(status["sibling_panels_scored"], ["plant_pt"])
        self.assertEqual(coverage["coverage_basis"],
                         cg.BASIS_TARGET_ASSIGNED)
        # Soft-fail still emits an empty gated FASTQ for channel typing.
        self.assertEqual((self.tmp / "gated.fastq.gz").read_bytes(), b"")

    def test_main_records_fallback_basis(self):
        """No bundle → basis is total_recruited, sample still gated."""
        def fake_run(cmd, stdout=None, **kwargs):
            self.assertEqual(cmd[0], "seqkit")
            return subprocess.CompletedProcess(
                cmd, 0,
                stdout="file\tformat\ttype\tnum_seqs\tsum_len\n"
                       "in\tFASTQ\tDNA\t2\t20000000\n")

        rc, status, coverage = self._run_main(self.tmp, fake_run)

        self.assertEqual(rc, 0)
        self.assertEqual(status["status"], "ok")
        self.assertEqual(status["coverage_basis"], cg.BASIS_TOTAL_RECRUITED)
        self.assertEqual(status["estimated_cov"], 50.0)
        self.assertEqual(status["sibling_panels_scored"], [])
        self.assertNotIn("target_assigned_bases", status)

    def test_main_zero_recruited_bases(self):
        """Nothing recruited → no sibling fraction to report, still exits 0."""
        def fake_run(cmd, stdout=None, **kwargs):
            return subprocess.CompletedProcess(
                cmd, 0,
                stdout="file\tformat\ttype\tnum_seqs\tsum_len\n"
                       "in\tFASTQ\tDNA\t0\t0\n")

        rc, status, _ = self._run_main(self.tmp, fake_run)

        self.assertEqual(rc, 0)
        self.assertEqual(status["status"], "low_coverage")
        self.assertEqual(status["estimated_cov"], 0.0)
        self.assertNotIn("sibling_organelle_fraction", status)

    def test_main_subsample_scales_target_estimate(self):
        """Over-MAX target depth subsamples; post cov reflects the target."""
        for panel in ("plant_mt", "plant_pt"):
            (self.tmp / f"{panel}.mmi").write_bytes(b"")
        gated = self.tmp / "gated.fastq.gz"

        def fake_run(cmd, stdout=None, **kwargs):
            if cmd[0] == "seqkit":
                # Second call measures the subsampled output.
                total = (10_000_000 if Path(cmd[3]).name.startswith("gated")
                         else 200_000_000)
                return subprocess.CompletedProcess(
                    cmd, 0,
                    stdout="file\tformat\ttype\tnum_seqs\tsum_len\n"
                           f"in\tFASTQ\tDNA\t2\t{total}\n")
            if cmd[0] == "bash":
                gated.write_bytes(b"")
                return subprocess.CompletedProcess(cmd, 0)
            panel = Path(cmd[-2]).stem
            if panel == "plant_mt":
                stdout.write(paf_line("mito", 160000000, 0, 150000000))
            else:
                stdout.write(paf_line("plastid", 40000000, 0, 39000000))
            return subprocess.CompletedProcess(cmd, 0)

        rc, status, coverage = self._run_main(self.tmp, fake_run)

        self.assertEqual(rc, 0)
        self.assertEqual(status["status"], "ok")
        # 160 Mb mito / 400 kb = 400x, over the 300x MAX.
        self.assertEqual(status["estimated_cov"], 400.0)
        self.assertTrue(coverage["subsampled"])
        self.assertEqual(coverage["fraction"], 0.75)
        # Realised reduction 10/200 = 0.05 → 400 × 0.05 = 20x.
        self.assertEqual(coverage["post_subsample_cov"], 20.0)


if __name__ == "__main__":
    unittest.main()
