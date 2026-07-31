"""Unit tests for bin/bin_target.py — spec §2 stage 10 decision branches.

Tests use synthetic in-memory data; no Nextflow or external tools required.
All coverage targets the pure select_primary() function plus the supporting
signal functions (longest_orf_aa, check_circularity).

Cases:
  1. Single obvious winner + 3 off-target contigs.
  2. Two high-cov contigs, only one passes ORF.
  3. Coverage tie broken by identity rank.
  4. No contig clears identity threshold → empty primaries.
  5. animal_mt synthetic circular contig → circular=True.
  6. animal_mt synthetic linear contig → circular=False.
  7. plant_mt multi-contig: all 3 candidates selected.
  8. animal_mt genetic-code trial: ORF only valid under table 5.
"""

import importlib.util
import unittest
from pathlib import Path

BIN_TARGET = Path(__file__).resolve().parents[2] / "bin" / "bin_target.py"

_spec = importlib.util.spec_from_file_location("bin_target", BIN_TARGET)
bt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bt)


# --- helpers ---

def _row(
    contig_id="ctg1",
    length_bp=17000,
    coverage=100.0,
    ref_identity_pct=95.0,
    aligned_frac=0.90,
    orf_ok=True,
    orf_aa_len=500,
    selected_genetic_code=2,
    classification="target_candidate",
):
    return {
        "contig_id": contig_id,
        "length_bp": length_bp,
        "coverage": coverage,
        "ref_identity_pct": ref_identity_pct,
        "aligned_frac": aligned_frac,
        "orf_ok": orf_ok,
        "orf_aa_len": orf_aa_len,
        "selected_genetic_code": selected_genetic_code,
        "classification": classification,
    }


def _off_target(contig_id, coverage=4.0):
    return _row(
        contig_id=contig_id,
        coverage=coverage,
        ref_identity_pct=10.0,
        aligned_frac=0.05,
        orf_ok=False,
        orf_aa_len=20,
        classification="off_target",
    )


# --- synthetic reference for mappy alignment tests ---

# Simple pseudo-random 1500 bp reference (deterministic, no external file).
_REF_UNIT = "ATGCTAGCTAGCTAGCATGCATGCATGCTAGCTAGCATGCATGCATGC"
REF_SEQ = (_REF_UNIT * 32)[:1500]
# High-identity query ≈ reference.
HIGH_ID_SEQ = REF_SEQ
# Off-target query: poly-N pattern, should not align to REF_SEQ.
OFF_TARGET_SEQ = "AAAAAAAAAA" * 150


class TestSelectPrimary(unittest.TestCase):

    def test_case1_single_winner(self):
        """1 target_candidate + 3 off_target → primary is the candidate."""
        rows = [
            _row("ctg_mt", classification="target_candidate"),
            _off_target("ctg_nuc1"),
            _off_target("ctg_nuc2"),
            _off_target("ctg_nuc3"),
        ]
        primaries, secondaries = bt.select_primary(rows, "animal_mt")
        self.assertEqual(len(primaries), 1)
        self.assertEqual(primaries[0]["contig_id"], "ctg_mt")
        self.assertEqual(len(secondaries), 3)
        self.assertTrue(all(r["classification"] == "off_target" for r in secondaries))

    def test_case2_orf_breaks_tie(self):
        """Two high-cov contigs; only one passes ORF → ORF-passer is primary."""
        rows = [
            _row("ctg_A", coverage=80.0, ref_identity_pct=90.0,
                 orf_ok=True, orf_aa_len=400, classification="target_candidate"),
            _row("ctg_B", coverage=80.0, ref_identity_pct=90.0,
                 orf_ok=False, orf_aa_len=50, classification="low_confidence"),
        ]
        primaries, secondaries = bt.select_primary(rows, "animal_mt")
        self.assertEqual(len(primaries), 1)
        self.assertEqual(primaries[0]["contig_id"], "ctg_A")
        self.assertEqual(len(secondaries), 1)
        self.assertEqual(secondaries[0]["contig_id"], "ctg_B")
        self.assertEqual(secondaries[0]["classification"], "low_confidence")

    def test_case3_coverage_tie_identity_breaks(self):
        """Coverage tie: higher identity wins."""
        rows = [
            _row("ctg_hi_id", coverage=100.0, ref_identity_pct=98.0,
                 orf_aa_len=400, classification="target_candidate"),
            _row("ctg_lo_id", coverage=100.0, ref_identity_pct=82.0,
                 orf_aa_len=400, classification="target_candidate"),
        ]
        primaries, secondaries = bt.select_primary(rows, "animal_mt")
        self.assertEqual(primaries[0]["contig_id"], "ctg_hi_id")
        sec_ids = [s["contig_id"] for s in secondaries]
        self.assertIn("ctg_lo_id", sec_ids)
        # Runner-up reclassified.
        reclassified = next(s for s in secondaries if s["contig_id"] == "ctg_lo_id")
        self.assertEqual(reclassified["classification"], "secondary_target")

    def test_case4_no_candidates(self):
        """No contig clears identity → empty primaries, all in secondaries."""
        rows = [
            _off_target("ctg1"),
            _off_target("ctg2"),
        ]
        primaries, secondaries = bt.select_primary(rows, "animal_mt")
        self.assertEqual(primaries, [])
        self.assertEqual(len(secondaries), 2)
        self.assertTrue(all(r["classification"] == "off_target" for r in secondaries))

    def test_case7_plant_mt_multi_contig(self):
        """plant_mt: all 3 target_candidates selected as primaries."""
        rows = [
            _row("ctg_mt1", coverage=60.0, classification="target_candidate"),
            _row("ctg_mt2", coverage=55.0, classification="target_candidate"),
            _row("ctg_mt3", coverage=50.0, classification="target_candidate"),
        ]
        primaries, secondaries = bt.select_primary(rows, "plant_mt")
        self.assertEqual(len(primaries), 3)
        self.assertEqual(secondaries, [])


class TestLongestOrfAa(unittest.TestCase):

    def _build_orf_seq(self, aa_count: int, table: int) -> str:
        """Build a nucleotide sequence with a known ORF of aa_count amino acids."""
        from Bio.Data.CodonTable import unambiguous_dna_by_id
        ct = unambiguous_dna_by_id[table]
        # Pick the first non-stop sense codon from the table.
        sense_codon = next(
            c for c, a in ct.forward_table.items()
            if a != "*" and c not in ct.stop_codons
        )
        stop_codon = ct.stop_codons[0]
        return "ATG" + (sense_codon * (aa_count - 1)) + stop_codon

    def test_basic_orf_detection(self):
        """Sequence with a 200 aa ORF is detected and length returned."""
        seq = self._build_orf_seq(200, table=1)
        length, selected = bt.longest_orf_aa(seq, [1])
        self.assertGreaterEqual(length, 200)

    def test_case8_orf_longer_under_table5(self):
        """Table 5 wins when the RC frame 0 is longer under t5 than t2.

        Sequence: ATG + (TCT)*300 + TAA
        The RC is: TTA + (AGA)*300 + CAT

        RC frame 0 under table 2: AGA=stop after first codon (TTA) → 1 aa run
        RC frame 0 under table 5: AGA=Ser → 302 aa uninterrupted run
        Forward frame 0 + other frames: ~301 aa (identical under t2 and t5
        since TCT=Ser in both).

        So best_t2 ≈ 301 aa; best_t5 = 302 aa (RC frame 0 adds one more).
        The trial [2, 5] selects table 5.
        """
        # TCT = Ser in tables 2 and 5; RC of TCT = AGA = stop(t2) / Ser(t5).
        seq = "ATG" + ("TCT" * 300) + "TAA"

        length_t2, _ = bt.longest_orf_aa(seq, [2])
        length_t5, table_t5 = bt.longest_orf_aa(seq, [5])

        # Table 5 produces a strictly longer ORF than table 2.
        self.assertGreater(length_t5, length_t2)
        self.assertEqual(table_t5, 5)

        # Combined trial selects table 5 (longer ORF wins).
        length_both, table_both = bt.longest_orf_aa(seq, [2, 5])
        self.assertEqual(table_both, 5)
        self.assertGreaterEqual(length_both, bt.MIN_ORF_AA)


class TestCheckCircularity(unittest.TestCase):

    def _circular_contig(self, length: int = 4000) -> str:
        """Build a contig whose prefix == suffix (simulates Flye circular overlap)."""
        unit = "ATGCTAGCTAGCGATCGATCGATCGATCTAGCATGCTAGCGATCGAT"
        overlap = (unit * 7)[:bt.END_OVERLAP_WINDOW]
        middle_len = length - 2 * len(overlap)
        middle = ("GCGCGCTAGCTAGCTAGCGCGCGCATGCATGCATGCTAGCTAGCTAGC" * 100)[:middle_len]
        return overlap + middle + overlap

    def _linear_contig(self, length: int = 4000) -> str:
        """Build a contig with completely different prefix and suffix."""
        prefix = ("ATGCTAGCTAGCGATCGATCGATCGATCTAGCATGCTAGCGATCGAT" * 7)[:bt.END_OVERLAP_WINDOW]
        suffix = ("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" * 7)[:bt.END_OVERLAP_WINDOW]
        middle_len = length - len(prefix) - len(suffix)
        middle = ("GCGCATGCTAGCGATCGCGCATGCATGCTAGCTAGCGCGCATGCATGC" * 100)[:middle_len]
        return prefix + middle + suffix

    def test_case5_circular_detected(self):
        """Contig with end-overlap detected as circular."""
        seq = self._circular_contig()
        self.assertTrue(bt.check_circularity(seq))

    def test_case6_linear_not_circular(self):
        """Contig without end-overlap detected as linear."""
        seq = self._linear_contig()
        self.assertFalse(bt.check_circularity(seq))

    def test_too_short_for_circularity(self):
        """Sequence shorter than 2 × END_OVERLAP_WINDOW returns False."""
        seq = "ATGC" * 10
        self.assertFalse(bt.check_circularity(seq))


if __name__ == "__main__":
    unittest.main()
