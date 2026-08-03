"""Unit tests for bin/bin_target.py — spec §2 stage 10, §3.7.

Tests use synthetic in-memory data; no Nextflow or external tools
required. Coverage targets the selection path: merge_intervals(),
align_to_ref(), classify_contigs(), select_primary(), load_panels(),
resolve_circularity() and the supporting signal functions.

Cases:
  1. Single obvious winner + 3 off-target contigs.
  2. Two contigs, only one a candidate → candidate is primary.
  3. Coverage tie broken by identity rank.
  4. No contig clears the homology floors → empty primaries.
  5. Synthetic end-overlap contig → check_circularity True.
  6. Synthetic linear contig → check_circularity False.
  7. plant_mt multi-contig: all 3 candidates selected.
  8. animal_mt genetic-code trial: ORF only valid under table 5.
  9. Interval merge: overlapping/adjacent/nested/disjoint/single/empty.
 10. Contig scoring higher against the sibling panel → sibling_organelle.
 11. Contig scoring higher against the declared panel → candidate.
 12. Sibling .mmi missing from the bundle → warns, still selects, 0.
 13. Flye circ. = Y, no end overlap → circular via flye_circ.
 14. Flye circ. = N, synthetic end overlap → circular via end_overlap.
 15. Flye circ. = N, no overlap → not circular.
 16. plant_mt with 2 mt + 2 plastid contigs → exactly the 2 mt emitted.
 17. All contigs below min_aligned_frac → empty target, rows classified.
 18. Per-target thresholds: admitted under plant_mt, rejected animal_mt.
 19. plant_mt low-coverage candidate flagged, not filtered (§5.2).
 20. plant_pt canonical graph + C3 selection → C4 path1 substituted.
 21. plant_pt canonical graph, C3 selected 0 → substitution withheld,
     empty target.fasta, isoforms retained (task 24 §4 case 2).
 22. Sibling-organelle carry-over fraction + warning (task 25 §3.1).
 22. plant_pt resolved_circle / non_canonical → no substitution.
"""

import importlib.util
import io
import json
import random
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
BIN_TARGET = BIN_DIR / "bin_target.py"

# bin_target.py imports plastid_canonicalise (C4) lazily on the
# plant_pt branch, relying on Nextflow's bin/ staging at runtime.
sys.path.insert(0, str(BIN_DIR))

_spec = importlib.util.spec_from_file_location("bin_target", BIN_TARGET)
bt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bt)

import mappy  # noqa: E402  (imported after the module under test loads)


# --- synthetic sequence helpers ---

def _random_dna(length: int, seed: int) -> str:
    """Deterministic pseudo-random DNA — complex enough for minimap2."""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


# Two unrelated 6 kb "organelle panels". A contig drawn from one does
# not align to the other, which is what the sibling test exploits.
REF_MT = _random_dna(6000, seed=11)
REF_PT = _random_dna(6000, seed=22)


def _mutate(seq: str, rate: float, seed: int) -> str:
    """Substitute `rate` of the bases — a diverged, still-alignable copy."""
    rng = random.Random(seed)
    return "".join(
        rng.choice([b for b in "ACGT" if b != base])
        if rng.random() < rate else base
        for base in seq
    )


def _canonical_gfa() -> str:
    """3-edge plastid graph (LSC/IR/SSC) that C4 calls `canonical`.

    Edge sequences are unrelated to REF_PT — C4 keys on graph topology
    alone, which is exactly the defect task 24 guards against.
    """
    edges = [
        ("edge_A", _random_dna(900, seed=901), 1.0),
        ("edge_B", _random_dna(250, seed=902), 2.0),
        ("edge_C", _random_dna(150, seed=903), 1.0),
    ]
    lines = ["H\tVN:Z:1.0"] + [
        f"S\t{name}\t{seq}\tdp:f:{depth}" for name, seq, depth in edges
    ]
    return "\n".join(lines) + "\n"


def _resolved_circle_gfa() -> str:
    """Single-edge graph — C4's `resolved_circle` branch."""
    return (
        "H\tVN:Z:1.0\n"
        f"S\tedge_A\t{_random_dna(1200, seed=904)}\tdp:f:1.0\n"
    )


def _mt_aligner():
    return mappy.Aligner(seq=REF_MT, preset="map-ont", best_n=5)


def _pt_aligner():
    return mappy.Aligner(seq=REF_PT, preset="map-ont", best_n=5)


def _thresholds(
    min_identity=80.0,
    min_aligned_frac=0.30,
    emit="single",
    max_contigs=1,
    low_coverage_fraction=0.05,
):
    return {
        "min_identity": min_identity,
        "min_aligned_frac": min_aligned_frac,
        "emit": emit,
        "max_contigs": max_contigs,
        "low_coverage_fraction": low_coverage_fraction,
    }


PLANT_MT_THRESHOLDS = _thresholds(
    min_identity=70.0, min_aligned_frac=0.15, emit="all", max_contigs=20)


def _info_entry(seq: str, coverage: float, flye_circular=False):
    return {
        "length": len(seq),
        "coverage": coverage,
        "flye_circular": flye_circular,
    }


# --- pre-classified row helper for the pure select_primary tests ---

def _row(
    contig_id="ctg1",
    length_bp=17000,
    coverage=100.0,
    ref_identity_pct=95.0,
    aligned_frac=0.90,
    orf_aa_len=500,
    classification="target_candidate",
):
    return {
        "contig_id": contig_id,
        "length_bp": length_bp,
        "coverage": coverage,
        "ref_identity_pct": ref_identity_pct,
        "aligned_frac": aligned_frac,
        "panel_aligned_frac": {"animal_mt": aligned_frac},
        "panel_identity_pct": {"animal_mt": ref_identity_pct},
        "winning_panel": "animal_mt",
        "circular": False,
        "circular_method": "none",
        "orf_ok": orf_aa_len >= bt.MIN_ORF_AA,
        "orf_aa_len": orf_aa_len,
        "selected_genetic_code": 2,
        "low_coverage_candidate": False,
        "classification": classification,
    }


def _off_target(contig_id, coverage=4.0):
    return _row(
        contig_id=contig_id,
        coverage=coverage,
        ref_identity_pct=10.0,
        aligned_frac=0.05,
        orf_aa_len=20,
        classification="off_target",
    )


class TestSelectPrimary(unittest.TestCase):

    def test_case1_single_winner(self):
        """1 target_candidate + 3 off_target → primary is the candidate."""
        rows = [
            _row("ctg_mt", classification="target_candidate"),
            _off_target("ctg_nuc1"),
            _off_target("ctg_nuc2"),
            _off_target("ctg_nuc3"),
        ]
        primaries, secondaries = bt.select_primary(
            rows, "animal_mt", _thresholds())
        self.assertEqual(len(primaries), 1)
        self.assertEqual(primaries[0]["contig_id"], "ctg_mt")
        self.assertEqual(len(secondaries), 3)
        self.assertTrue(
            all(r["classification"] == "off_target" for r in secondaries))

    def test_case2_non_candidate_not_selected(self):
        """Equal coverage; only the candidate row is eligible.

        Re-baselined for §3.7.5 — ORF integrity no longer gates
        selection, so the loser here is off_target on homology.
        """
        rows = [
            _row("ctg_A", coverage=80.0, ref_identity_pct=90.0,
                 classification="target_candidate"),
            _row("ctg_B", coverage=80.0, ref_identity_pct=90.0,
                 aligned_frac=0.05, classification="off_target"),
        ]
        primaries, secondaries = bt.select_primary(
            rows, "animal_mt", _thresholds())
        self.assertEqual(len(primaries), 1)
        self.assertEqual(primaries[0]["contig_id"], "ctg_A")
        self.assertEqual(len(secondaries), 1)
        self.assertEqual(secondaries[0]["contig_id"], "ctg_B")
        self.assertEqual(secondaries[0]["classification"], "off_target")

    def test_case3_coverage_tie_identity_breaks(self):
        """Coverage tie: higher identity wins (no coverage-spike gate)."""
        rows = [
            _row("ctg_hi_id", coverage=100.0, ref_identity_pct=98.0),
            _row("ctg_lo_id", coverage=100.0, ref_identity_pct=82.0),
        ]
        primaries, secondaries = bt.select_primary(
            rows, "animal_mt", _thresholds())
        self.assertEqual(primaries[0]["contig_id"], "ctg_hi_id")
        reclassified = next(
            s for s in secondaries if s["contig_id"] == "ctg_lo_id")
        self.assertEqual(reclassified["classification"], "secondary_target")

    def test_case3b_identity_tie_frac_breaks(self):
        """Coverage and identity tie: higher aligned fraction wins."""
        rows = [
            _row("ctg_hi_frac", ref_identity_pct=90.0, aligned_frac=0.95),
            _row("ctg_lo_frac", ref_identity_pct=90.0, aligned_frac=0.40),
        ]
        primaries, _ = bt.select_primary(rows, "animal_mt", _thresholds())
        self.assertEqual(primaries[0]["contig_id"], "ctg_hi_frac")

    def test_case4_no_candidates(self):
        """No contig clears the floors → empty primaries, all secondary."""
        rows = [_off_target("ctg1"), _off_target("ctg2")]
        primaries, secondaries = bt.select_primary(
            rows, "animal_mt", _thresholds())
        self.assertEqual(primaries, [])
        self.assertEqual(len(secondaries), 2)

    def test_case7_plant_mt_multi_contig(self):
        """plant_mt: all 3 target_candidates selected as primaries."""
        rows = [
            _row("ctg_mt1", coverage=60.0),
            _row("ctg_mt2", coverage=55.0),
            _row("ctg_mt3", coverage=50.0),
        ]
        primaries, secondaries = bt.select_primary(
            rows, "plant_mt", PLANT_MT_THRESHOLDS)
        self.assertEqual(len(primaries), 3)
        self.assertEqual(secondaries, [])

    def test_max_contigs_cap(self):
        """emit=all honours max_contigs; the surplus becomes secondary."""
        rows = [
            _row(f"ctg_mt{i}", coverage=float(60 - i)) for i in range(5)
        ]
        thresholds = _thresholds(
            emit="all", max_contigs=2, min_aligned_frac=0.15)
        primaries, secondaries = bt.select_primary(
            rows, "plant_mt", thresholds)
        self.assertEqual(
            [r["contig_id"] for r in primaries], ["ctg_mt0", "ctg_mt1"])
        self.assertEqual(len(secondaries), 3)
        self.assertTrue(all(
            r["classification"] == "secondary_target"
            for r in secondaries))

    def test_case19_low_coverage_candidate_flagged(self):
        """§5.2: a 0.01× candidate is emitted, but flagged."""
        rows = [
            _row("ctg_mt_top", coverage=150.0),
            _row("ctg_mt_numt", coverage=1.5),
        ]
        primaries, _ = bt.select_primary(
            rows, "plant_mt", PLANT_MT_THRESHOLDS)
        flags = {r["contig_id"]: r["low_coverage_candidate"]
                 for r in primaries}
        self.assertEqual(len(primaries), 2)
        self.assertFalse(flags["ctg_mt_top"])
        self.assertTrue(flags["ctg_mt_numt"])


class TestMergeIntervals(unittest.TestCase):
    """Case 9 — the merged-interval metric of spec §3.7.2."""

    def _merged_len(self, intervals):
        return sum(e - s for s, e in bt.merge_intervals(intervals))

    def test_empty(self):
        self.assertEqual(bt.merge_intervals([]), [])

    def test_single(self):
        self.assertEqual(bt.merge_intervals([(10, 50)]), [(10, 50)])

    def test_overlapping(self):
        self.assertEqual(
            bt.merge_intervals([(0, 100), (50, 150)]), [(0, 150)])

    def test_adjacent(self):
        self.assertEqual(
            bt.merge_intervals([(0, 100), (100, 200)]), [(0, 200)])

    def test_nested(self):
        self.assertEqual(
            bt.merge_intervals([(0, 500), (100, 200)]), [(0, 500)])

    def test_disjoint(self):
        self.assertEqual(
            bt.merge_intervals([(300, 400), (0, 100)]),
            [(0, 100), (300, 400)])
        self.assertEqual(self._merged_len([(300, 400), (0, 100)]), 200)

    def test_unsorted_mixed(self):
        """Real minimap2 output: unsorted, partly overlapping blocks."""
        intervals = [(900, 1000), (0, 200), (150, 400), (400, 500)]
        self.assertEqual(
            bt.merge_intervals(intervals), [(0, 500), (900, 1000)])
        self.assertEqual(self._merged_len(intervals), 600)


class TestAlignToRef(unittest.TestCase):

    def test_no_hits(self):
        """Unrelated sequence → (0.0, 0.0)."""
        identity, frac = bt.align_to_ref(REF_PT[:2000], _mt_aligner())
        self.assertEqual((identity, frac), (0.0, 0.0))

    def test_empty_sequence(self):
        self.assertEqual(bt.align_to_ref("", _mt_aligner()), (0.0, 0.0))

    def test_full_length_hit(self):
        """Self-alignment → near-total merged fraction, high identity."""
        identity, frac = bt.align_to_ref(REF_MT[:4000], _mt_aligner())
        self.assertGreater(frac, 0.95)
        self.assertGreater(identity, 95.0)

    def test_merged_exceeds_single_block(self):
        """Two separated hits merge to roughly the sum of their lengths."""
        seq = REF_MT[:1500] + _random_dna(500, seed=99) + REF_MT[4000:5500]
        _, frac = bt.align_to_ref(seq, _mt_aligner())
        self.assertGreater(frac, 0.75)


class TestClassifyContigs(unittest.TestCase):

    def _classify(self, sequences, target, aligners, thresholds,
                  coverages=None, circ=None):
        coverages = coverages or {}
        circ = circ or {}
        info = {
            cid: _info_entry(
                seq, coverages.get(cid, 50.0), circ.get(cid, False))
            for cid, seq in sequences.items()
        }
        return bt.classify_contigs(
            info, sequences, aligners, target, [1], thresholds)

    def test_case10_sibling_organelle(self):
        """Contig scoring higher against the sibling panel is excluded."""
        sequences = {"ctg_plastid": REF_PT[:4000]}
        rows = self._classify(
            sequences, "plant_mt",
            {"plant_mt": _mt_aligner(), "plant_pt": _pt_aligner()},
            PLANT_MT_THRESHOLDS,
        )
        row = rows[0]
        self.assertEqual(row["classification"], "sibling_organelle")
        self.assertEqual(row["winning_panel"], "plant_pt")
        self.assertGreater(row["panel_aligned_frac"]["plant_pt"], 0.9)
        primaries, _ = bt.select_primary(
            rows, "plant_mt", PLANT_MT_THRESHOLDS)
        self.assertEqual(primaries, [])

    def test_case11_declared_panel_wins(self):
        """Both panels non-zero, declared higher → target_candidate."""
        # 3 kb of mt sequence + 1 kb of plastid sequence: the declared
        # panel wins on merged fraction (~0.75 vs ~0.25).
        sequences = {"ctg_mt": REF_MT[:3000] + REF_PT[:1000]}
        rows = self._classify(
            sequences, "plant_mt",
            {"plant_mt": _mt_aligner(), "plant_pt": _pt_aligner()},
            PLANT_MT_THRESHOLDS,
        )
        row = rows[0]
        self.assertGreater(row["panel_aligned_frac"]["plant_pt"], 0.0)
        self.assertGreater(
            row["panel_aligned_frac"]["plant_mt"],
            row["panel_aligned_frac"]["plant_pt"])
        self.assertEqual(row["classification"], "target_candidate")
        self.assertEqual(row["winning_panel"], "plant_mt")

    def test_case16_plant_mt_two_mt_two_plastid(self):
        """Regression test for the defect motivating task 23.

        The plastid contigs carry the *highest* coverage, exactly as in
        the INT-PLANT-01-mt fixture. Only the mt contigs may be emitted.
        """
        sequences = {
            "ctg_mt_1": REF_MT[:3000],
            "ctg_mt_2": REF_MT[3000:6000],
            "ctg_pt_1": REF_PT[:3000],
            "ctg_pt_2": REF_PT[3000:6000],
        }
        coverages = {
            "ctg_mt_1": 7.0, "ctg_mt_2": 5.0,
            "ctg_pt_1": 168.0, "ctg_pt_2": 137.0,
        }
        rows = self._classify(
            sequences, "plant_mt",
            {"plant_mt": _mt_aligner(), "plant_pt": _pt_aligner()},
            PLANT_MT_THRESHOLDS, coverages=coverages,
        )
        primaries, secondaries = bt.select_primary(
            rows, "plant_mt", PLANT_MT_THRESHOLDS)
        self.assertEqual(
            sorted(r["contig_id"] for r in primaries),
            ["ctg_mt_1", "ctg_mt_2"])
        self.assertTrue(all(
            r["classification"] == "sibling_organelle"
            for r in secondaries))

    def test_case17_all_below_min_aligned_frac(self):
        """Nothing clears the floor → no primaries, every row classified."""
        sequences = {
            "ctg1": _random_dna(3000, seed=101),
            "ctg2": _random_dna(3000, seed=102),
        }
        rows = self._classify(
            sequences, "animal_mt", {"animal_mt": _mt_aligner()},
            _thresholds(),
        )
        self.assertTrue(
            all(r["classification"] == "off_target" for r in rows))
        self.assertTrue(all(r["winning_panel"] == "none" for r in rows))
        primaries, secondaries = bt.select_primary(
            rows, "animal_mt", _thresholds())
        self.assertEqual(primaries, [])
        self.assertEqual(len(secondaries), 2)

    def test_case18_per_target_thresholds(self):
        """Same contig: admitted under plant_mt bounds, rejected animal.

        1 kb of panel homology in a 4 kb contig → merged fraction ~0.25,
        which clears plant_mt's 0.15 floor but not animal_mt's 0.30.
        """
        sequences = {"ctg": REF_MT[:1000] + _random_dna(3000, seed=77)}
        aligners = {"plant_mt": _mt_aligner()}
        rows = bt.classify_contigs(
            {"ctg": _info_entry(sequences["ctg"], 50.0)},
            sequences, aligners, "plant_mt", [1], PLANT_MT_THRESHOLDS)
        self.assertGreater(rows[0]["aligned_frac"], 0.15)
        self.assertLess(rows[0]["aligned_frac"], 0.30)
        self.assertEqual(rows[0]["classification"], "target_candidate")

        rows = bt.classify_contigs(
            {"ctg": _info_entry(sequences["ctg"], 50.0)},
            sequences, {"animal_mt": _mt_aligner()}, "animal_mt", [2],
            _thresholds())
        self.assertEqual(rows[0]["classification"], "off_target")

    def test_identity_floor_rejects(self):
        """Aligned fraction clears the floor but identity does not."""
        thresholds = _thresholds(min_identity=99.0, min_aligned_frac=0.30)
        sequences = {"ctg": _mutate(REF_MT[:3000], rate=0.05, seed=55)}
        rows = self._classify(
            sequences, "animal_mt", {"animal_mt": _mt_aligner()},
            thresholds)
        self.assertGreater(rows[0]["aligned_frac"], 0.30)
        self.assertEqual(rows[0]["classification"], "off_target")


class TestCircularity(unittest.TestCase):

    def _circular_contig(self, length: int = 4000) -> str:
        """Contig whose prefix == suffix (un-trimmed Flye overlap)."""
        overlap = _random_dna(bt.END_OVERLAP_WINDOW, seed=31)
        middle = _random_dna(length - 2 * len(overlap), seed=32)
        return overlap + middle + overlap

    def _linear_contig(self, length: int = 4000) -> str:
        return _random_dna(length, seed=41)

    def test_case5_circular_detected(self):
        self.assertTrue(bt.check_circularity(self._circular_contig()))

    def test_case6_linear_not_circular(self):
        self.assertFalse(bt.check_circularity(self._linear_contig()))

    def test_too_short_for_circularity(self):
        """Sequence shorter than 2 × END_OVERLAP_WINDOW returns False."""
        self.assertFalse(bt.check_circularity("ATGC" * 10))

    def test_case13_flye_circ_yes(self):
        """circ. = Y with no end overlap → circular via flye_circ."""
        seq = self._linear_contig()
        self.assertEqual(
            bt.resolve_circularity(True, seq), (True, "flye_circ"))

    def test_case14_end_overlap_fallback(self):
        """circ. = N with a 300 bp end overlap → end_overlap."""
        seq = self._circular_contig()
        self.assertEqual(
            bt.resolve_circularity(False, seq), (True, "end_overlap"))

    def test_case15_not_circular(self):
        """circ. = N with no overlap → none."""
        seq = self._linear_contig()
        self.assertEqual(
            bt.resolve_circularity(False, seq), (False, "none"))

    def test_circularity_recorded_per_contig(self):
        """classify_contigs records both fields for every contig."""
        seq = self._linear_contig()
        rows = bt.classify_contigs(
            {"ctg": _info_entry(seq, 50.0, flye_circular=True)},
            {"ctg": seq}, {"animal_mt": _mt_aligner()}, "animal_mt",
            [2], _thresholds())
        self.assertTrue(rows[0]["circular"])
        self.assertEqual(rows[0]["circular_method"], "flye_circ")


class TestParseAssemblyInfo(unittest.TestCase):

    def _write(self, body: str) -> Path:
        path = Path(self.tmp.name) / "assembly_info.txt"
        path.write_text(body)
        return path

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_reads_circ_column(self):
        path = self._write(
            "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\n"
            "contig_8\t16952\t71\tY\tN\t1\n"
            "contig_9\t16942\t31\tN\tN\t1\n"
        )
        info = bt.parse_assembly_info(path)
        self.assertTrue(info["contig_8"]["flye_circular"])
        self.assertFalse(info["contig_9"]["flye_circular"])
        self.assertEqual(info["contig_8"]["length"], 16952)
        self.assertEqual(info["contig_8"]["coverage"], 71.0)

    def test_tolerates_missing_circ_column(self):
        """Truncated rows still parse; short rows are skipped."""
        path = self._write(
            "#seq_name\tlength\tcov.\n"
            "contig_1\t1000\t10\n"
            "junk\t5\n"
        )
        info = bt.parse_assembly_info(path)
        self.assertEqual(list(info), ["contig_1"])
        self.assertFalse(info["contig_1"]["flye_circular"])


class TestLoadPanels(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.ref_dir = Path(self.tmp.name)

    def _write_panel(self, name: str, seq: str) -> None:
        # mappy accepts FASTA in place of a prebuilt index; the .mmi
        # name is what the bundle layout uses.
        (self.ref_dir / f"{name}.mmi").write_text(f">{name}\n{seq}\n")

    def test_loads_declared_and_sibling(self):
        self._write_panel("plant_mt", REF_MT)
        self._write_panel("plant_pt", REF_PT)
        aligners, scored = bt.load_panels(self.ref_dir, "plant_mt")
        self.assertEqual(sorted(aligners), ["plant_mt", "plant_pt"])
        self.assertEqual(scored, ["plant_pt"])

    def test_animal_mt_has_no_sibling(self):
        self._write_panel("animal_mt", REF_MT)
        aligners, scored = bt.load_panels(self.ref_dir, "animal_mt")
        self.assertEqual(list(aligners), ["animal_mt"])
        self.assertEqual(scored, [])

    def test_case12_missing_sibling_warns_and_continues(self):
        """Older bundle without plant_pt.mmi → warn, do not fail."""
        self._write_panel("plant_mt", REF_MT)
        err = io.StringIO()
        with redirect_stderr(err):
            aligners, scored = bt.load_panels(self.ref_dir, "plant_mt")
        self.assertEqual(scored, [])
        self.assertIn("WARNING", err.getvalue())
        self.assertIn("plant_pt.mmi", err.getvalue())

        # Selection still works — the declared panel alone.
        seq = REF_MT[:3000]
        rows = bt.classify_contigs(
            {"ctg": _info_entry(seq, 50.0)}, {"ctg": seq}, aligners,
            "plant_mt", [1], PLANT_MT_THRESHOLDS)
        primaries, _ = bt.select_primary(
            rows, "plant_mt", PLANT_MT_THRESHOLDS)
        self.assertEqual(len(primaries), 1)

    def test_unloadable_sibling_warns(self):
        """A present but corrupt sibling index degrades gracefully."""
        self._write_panel("plant_mt", REF_MT)
        (self.ref_dir / "plant_pt.mmi").write_bytes(b"\x00not-an-index")
        err = io.StringIO()
        with redirect_stderr(err):
            aligners, scored = bt.load_panels(self.ref_dir, "plant_mt")
        self.assertEqual(scored, [])
        self.assertEqual(list(aligners), ["plant_mt"])
        self.assertIn("WARNING", err.getvalue())

    def test_missing_declared_panel_is_fatal(self):
        with self.assertRaises(SystemExit):
            bt.load_panels(self.ref_dir, "plant_mt")


class TestMainEndToEnd(unittest.TestCase):
    """Exercises the CLI wiring, outputs and metadata contract."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.refs = self.dir / "refs"
        self.refs.mkdir()

    def _run(
        self,
        target,
        sequences,
        coverages,
        circ=None,
        gfa_text=None,
        **overrides,
    ):
        circ = circ or {}
        (self.refs / "plant_mt.mmi").write_text(f">mt\n{REF_MT}\n")
        (self.refs / "plant_pt.mmi").write_text(f">pt\n{REF_PT}\n")
        (self.refs / "animal_mt.mmi").write_text(f">amt\n{REF_MT}\n")

        assembly = self.dir / "assembly.fasta"
        assembly.write_text("".join(
            f">{cid}\n{seq}\n" for cid, seq in sequences.items()))
        info = self.dir / "assembly_info.txt"
        info.write_text(
            "#seq_name\tlength\tcov.\tcirc.\n" + "".join(
                f"{cid}\t{len(seq)}\t{coverages[cid]}\t"
                f"{'Y' if circ.get(cid) else 'N'}\n"
                for cid, seq in sequences.items()))
        gfa = self.dir / "assembly.gfa"
        gfa.write_text(gfa_text or "H\tVN:Z:1.0\n")

        thresholds = _thresholds(
            min_identity=70.0, min_aligned_frac=0.15, emit="all",
            max_contigs=20)
        thresholds.update(overrides)

        argv = [
            "bin_target.py",
            "--assembly", str(assembly),
            "--assembly-info", str(info),
            "--gfa", str(gfa),
            "--ref-dir", str(self.refs),
            "--sample-id", "TEST-01",
            "--assembly-target", target,
            "--genetic-codes", "1",
            "--min-identity", str(thresholds["min_identity"]),
            "--min-aligned-frac", str(thresholds["min_aligned_frac"]),
            "--emit", thresholds["emit"],
            "--max-contigs", str(thresholds["max_contigs"]),
            "--low-coverage-fraction",
            str(thresholds["low_coverage_fraction"]),
            "--sibling-warn-fraction",
            str(thresholds.get("sibling_warn_fraction", 0.30)),
            "--out-target", str(self.dir / "target.fasta"),
            "--out-secondaries", str(self.dir / "secondaries.tsv"),
            "--out-metadata", str(self.dir / "bin_metadata.json"),
        ]
        original = sys.argv
        sys.argv = argv
        try:
            rc = bt.main()
        finally:
            sys.argv = original
        return (
            rc,
            (self.dir / "target.fasta").read_text(),
            (self.dir / "secondaries.tsv").read_text(),
            json.loads((self.dir / "bin_metadata.json").read_text()),
        )

    def test_plant_mt_excludes_plastid(self):
        """End-to-end §3.7.1: plastid contigs never reach target.fasta."""
        sequences = {
            "contig_3": REF_MT[:3000],
            "contig_5": REF_MT[3000:6000],
            "contig_1": REF_PT[:3000],
            "contig_6": REF_PT[3000:6000],
        }
        coverages = {
            "contig_3": 7.0, "contig_5": 5.0,
            "contig_1": 168.0, "contig_6": 137.0,
        }
        rc, target, secondaries, meta = self._run(
            "plant_mt", sequences, coverages)
        self.assertEqual(rc, 0)
        self.assertIn(">contig_3", target)
        self.assertIn(">contig_5", target)
        self.assertNotIn(">contig_1", target)
        self.assertNotIn(">contig_6", target)
        self.assertEqual(meta["n_target_selected"], 2)
        self.assertEqual(meta["sibling_panels_scored"], ["plant_pt"])
        self.assertEqual(
            meta["thresholds_applied"]["min_aligned_frac"], 0.15)

        # Every panel scored is reported per contig, in both outputs.
        header = secondaries.splitlines()[0].split("\t")
        self.assertIn("aligned_frac_plant_mt", header)
        self.assertIn("aligned_frac_plant_pt", header)
        audited = {c["contig_id"]: c for c in meta["contigs"]}
        self.assertEqual(len(audited), 4)
        self.assertEqual(
            audited["contig_1"]["classification"], "sibling_organelle")
        self.assertEqual(
            audited["contig_1"]["winning_panel"], "plant_pt")
        self.assertIn("plant_pt", audited["contig_3"]["panel_aligned_frac"])

    def test_animal_mt_circularity_and_single_emit(self):
        sequences = {
            "contig_8": REF_MT[:4000],
            "contig_9": REF_MT[2000:5000],
        }
        rc, target, _, meta = self._run(
            "animal_mt", sequences,
            {"contig_8": 71.0, "contig_9": 31.0},
            circ={"contig_8": True},
            emit="single", max_contigs=1,
        )
        self.assertEqual(rc, 0)
        self.assertEqual(target.count(">"), 1)
        self.assertIn(">contig_8", target)
        self.assertTrue(meta["circular"])
        self.assertEqual(meta["circular_method"], "flye_circ")
        self.assertEqual(
            meta["plastid_canonicalisation"]["branch"], "not_applicable")

    def _run_plant_pt(self, sequences, coverages, gfa_text=None):
        """plant_pt end-to-end with single-emit binning thresholds."""
        return self._run(
            "plant_pt", sequences, coverages, gfa_text=gfa_text,
            emit="single", max_contigs=1, min_identity=75.0,
            min_aligned_frac=0.30,
        )

    def test_plant_pt_selects_on_own_evidence(self):
        """Case 4 (task 24 §4): non_canonical graph, C3 selected ≥ 1.

        The stub GFA has no edges, so C4 reports non_canonical: C3's own
        selection is emitted unchanged and no isoforms are retained.
        """
        rc, target, _, meta = self._run_plant_pt(
            {"contig_2": REF_PT[:4000]}, {"contig_2": 148.0})
        self.assertEqual(rc, 0)
        self.assertIn(">contig_2", target)
        self.assertEqual(meta["n_target_selected"], 1)
        self.assertEqual(meta["sibling_panels_scored"], ["plant_mt"])
        self.assertEqual(meta["target_source"], bt.TARGET_SOURCE_C3)
        canon = meta["plastid_canonicalisation"]
        self.assertEqual(canon["branch"], "non_canonical")
        self.assertFalse(canon["substitution_applied"])
        self.assertNotIn("substitution_withheld_reason", canon)
        self.assertFalse((self.dir / "plastid_isoforms").exists())

    def test_plant_pt_canonical_with_selection_substitutes(self):
        """Case 1 (task 24 §4): canonical graph + C3 selection.

        The substitution proceeds and target.fasta becomes C4's path1,
        recorded as such in `target_source`.
        """
        rc, target, _, meta = self._run_plant_pt(
            {"contig_2": REF_PT[:4000]}, {"contig_2": 148.0},
            gfa_text=_canonical_gfa())
        self.assertEqual(rc, 0)
        self.assertEqual(meta["n_target_selected"], 1)
        self.assertEqual(meta["contigs_selected"], ["contig_2"])
        self.assertEqual(meta["target_source"], bt.TARGET_SOURCE_C4_PATH1)
        canon = meta["plastid_canonicalisation"]
        self.assertEqual(canon["branch"], "canonical")
        self.assertTrue(canon["substitution_applied"])
        self.assertNotIn("substitution_withheld_reason", canon)

        iso_dir = self.dir / "plastid_isoforms"
        path1 = iso_dir / "path1.fasta"
        self.assertTrue((iso_dir / "path2.fasta").exists())
        self.assertEqual(target, path1.read_text())
        self.assertNotIn(">contig_2", target)

    def test_plant_pt_canonical_without_selection_withholds(self):
        """Case 2 (task 24 §4) — the regression test for the defect.

        A 3-edge graph is a structural observation; with no contig
        passing target binning the substitution is withheld and the
        sample stays a visible negative.
        """
        rc, target, _, meta = self._run_plant_pt(
            {"contig_1": _random_dna(4000, seed=404)}, {"contig_1": 150.0},
            gfa_text=_canonical_gfa())
        self.assertEqual(rc, 0)
        self.assertEqual(target, "")
        self.assertEqual(meta["n_target_selected"], 0)
        self.assertEqual(meta["target_source"], bt.TARGET_SOURCE_C3)
        canon = meta["plastid_canonicalisation"]
        self.assertEqual(canon["branch"], "canonical")
        self.assertFalse(canon["substitution_applied"])
        self.assertEqual(
            canon["substitution_withheld_reason"],
            bt.SUBSTITUTION_WITHHELD_NO_SELECTION)

        # Isoforms stay on disk as operator diagnostics (task 24 §3).
        iso_dir = self.dir / "plastid_isoforms"
        self.assertTrue((iso_dir / "path1.fasta").exists())
        self.assertTrue((iso_dir / "path2.fasta").exists())

    def test_plant_pt_non_canonical_clears_stale_isoforms(self):
        """Off the canonical branch a pre-existing isoform dir is cleared.

        Only the withheld *canonical* case retains isoforms; a populated
        directory on any other branch is upstream residue.
        """
        stale = self.dir / "plastid_isoforms"
        stale.mkdir()
        (stale / "path1.fasta").write_text(">stale\nACGT\n")

        rc, _, _, meta = self._run_plant_pt(
            {"contig_2": REF_PT[:4000]}, {"contig_2": 148.0})
        self.assertEqual(rc, 0)
        self.assertEqual(
            meta["plastid_canonicalisation"]["branch"], "non_canonical")
        self.assertFalse(stale.exists())

    def test_plant_pt_resolved_circle_without_selection(self):
        """Case 3 (task 24 §4): resolved_circle, C3 selected 0.

        Unchanged behaviour — no substitution, no isoform directory.
        """
        rc, target, _, meta = self._run_plant_pt(
            {"contig_1": _random_dna(4000, seed=505)}, {"contig_1": 150.0},
            gfa_text=_resolved_circle_gfa())
        self.assertEqual(rc, 0)
        self.assertEqual(target, "")
        self.assertEqual(meta["n_target_selected"], 0)
        canon = meta["plastid_canonicalisation"]
        self.assertEqual(canon["branch"], "resolved_circle")
        self.assertFalse(canon["substitution_applied"])
        self.assertNotIn("substitution_withheld_reason", canon)
        self.assertFalse((self.dir / "plastid_isoforms").exists())

    def test_empty_selection_exits_zero(self):
        """Nothing selected → empty target.fasta, exit 0, rows recorded."""
        sequences = {"contig_1": _random_dna(3000, seed=303)}
        rc, target, secondaries, meta = self._run(
            "plant_mt", sequences, {"contig_1": 20.0})
        self.assertEqual(rc, 0)
        self.assertEqual(target, "")
        self.assertEqual(meta["n_target_selected"], 0)
        self.assertEqual(len(secondaries.splitlines()), 2)
        self.assertEqual(
            meta["contigs"][0]["classification"], "off_target")


class TestLongestOrfAa(unittest.TestCase):

    def _build_orf_seq(self, aa_count: int, table: int) -> str:
        """Build a nucleotide sequence with a known ORF of aa_count aa."""
        from Bio.Data.CodonTable import unambiguous_dna_by_id
        ct = unambiguous_dna_by_id[table]
        sense_codon = next(
            c for c, a in ct.forward_table.items()
            if a != "*" and c not in ct.stop_codons
        )
        stop_codon = ct.stop_codons[0]
        return "ATG" + (sense_codon * (aa_count - 1)) + stop_codon

    def test_basic_orf_detection(self):
        """Sequence with a 200 aa ORF is detected and length returned."""
        seq = self._build_orf_seq(200, table=1)
        length, _ = bt.longest_orf_aa(seq, [1])
        self.assertGreaterEqual(length, 200)

    def test_empty_sequence(self):
        """Zero-length input returns the fallback table, length 0."""
        self.assertEqual(bt.longest_orf_aa("", [1]), (0, 1))

    def test_case8_orf_longer_under_table5(self):
        """Table 5 wins when the RC frame 0 is longer under t5 than t2.

        Sequence: ATG + (TCT)*300 + TAA; RC: TTA + (AGA)*300 + CAT.
        AGA is a stop in table 2 but Ser in table 5, so the RC frame 0
        run is 1 aa under t2 and 302 aa under t5.
        """
        seq = "ATG" + ("TCT" * 300) + "TAA"

        length_t2, _ = bt.longest_orf_aa(seq, [2])
        length_t5, table_t5 = bt.longest_orf_aa(seq, [5])
        self.assertGreater(length_t5, length_t2)
        self.assertEqual(table_t5, 5)

        length_both, table_both = bt.longest_orf_aa(seq, [2, 5])
        self.assertEqual(table_both, 5)
        self.assertGreaterEqual(length_both, bt.MIN_ORF_AA)


class TestSiblingCarryover(unittest.TestCase):
    """Case 22 — task 25 §3.1, spec §2.1.5."""

    @staticmethod
    def _row(length, classification):
        return {"length_bp": length, "classification": classification}

    def test_fraction_and_warning_above_threshold(self):
        """The INT-PLANT-01-mt shape: most assembled bases are plastid."""
        rows = [
            self._row(86004, "sibling_organelle"),
            self._row(69287, "sibling_organelle"),
            self._row(68658, "target_candidate"),
            self._row(31797, "target_candidate"),
        ]
        summary = bt.sibling_organelle_summary(rows, 0.30)
        self.assertEqual(summary["assembly_bases"], 255746)
        self.assertEqual(summary["sibling_organelle_bases"], 155291)
        self.assertEqual(summary["sibling_organelle_fraction"], 0.6072)
        self.assertTrue(summary["sibling_carryover_warning"])

    def test_no_warning_below_threshold(self):
        rows = [
            self._row(10000, "sibling_organelle"),
            self._row(90000, "target_candidate"),
        ]
        summary = bt.sibling_organelle_summary(rows, 0.30)
        self.assertEqual(summary["sibling_organelle_fraction"], 0.1)
        self.assertFalse(summary["sibling_carryover_warning"])

    def test_no_siblings_scored(self):
        """animal_mt: nothing can bin as a sibling, so never warns."""
        summary = bt.sibling_organelle_summary(
            [self._row(16952, "target_candidate")], 0.30)
        self.assertEqual(summary["sibling_organelle_bases"], 0)
        self.assertEqual(summary["sibling_organelle_fraction"], 0.0)
        self.assertFalse(summary["sibling_carryover_warning"])

    def test_empty_assembly(self):
        """No contigs → no division by zero, no warning."""
        summary = bt.sibling_organelle_summary([], 0.30)
        self.assertEqual(summary["assembly_bases"], 0)
        self.assertEqual(summary["sibling_organelle_fraction"], 0.0)
        self.assertFalse(summary["sibling_carryover_warning"])


if __name__ == "__main__":
    unittest.main()
