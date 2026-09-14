"""Unit tests for scripts/refdata/build_proteins.py's rank-aware
breadth-first selection algorithm (task 48 §3, §6).

`scripts/refdata/` has no prior test coverage (§6) and is not `bin/`,
so `scripts/pytest.sh`'s coverage `--include` deliberately does not
report on it (task 48 §6) — these tests still run under that script,
just outside the coverage summary.

Cases, numbered to match task 48 §6:
  1. Breadth before depth.
  2. Recursion descends (class -> order/family -> genus).
  3. One record per genus at the leaf, longest translation wins.
  4. Under-fill.
  5. Unspent redistribution.
  6. Stratum disjointness.
  7. Determinism under shuffled input.
  8. Ragged lineages.
  9. `min_per_locus` interaction.
"""

import importlib.util
import random
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[2]
    / "scripts" / "refdata" / "build_proteins.py"
)
_spec = importlib.util.spec_from_file_location(
    "build_proteins", MODULE_PATH
)
bp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bp)


def _rec(accession, organism, lineage, gene="COX1", translation=None):
    return {
        "accession": accession,
        "organism": organism,
        "gene": gene,
        "translation": translation or ("M" * 100),
        "lineage": lineage,
    }


class TestBreadthBeforeDepth(unittest.TestCase):
    """§6.1 — 100 records from one insect order + 1 from Arachnida;
    budget 2 returns one of each, not two beetles."""

    def test_one_from_each_class(self):
        beetles = [
            _rec(f"BEETLE{i}", f"Genus{i} species", ["Insecta"])
            for i in range(100)
        ]
        mite = [_rec("MITE1", "Varroa destructor", ["Arachnida"])]
        items = [(r, bp._rel_lineage(r, None)) for r in beetles + mite]

        chosen = bp.select(items, 2, set(), set())

        self.assertEqual(len(chosen), 2)
        organisms = {r["organism"] for r in chosen}
        self.assertIn("Varroa destructor", organisms)
        beetle_count = sum(1 for r in chosen if r["accession"] != "MITE1")
        self.assertEqual(beetle_count, 1)


class TestRecursionDescends(unittest.TestCase):
    """§6.2 — budget large enough to exhaust order-level diversity
    spills into family, then genus."""

    def _build_pool(self):
        records = []
        # Coleoptera: 2 families, 5 genera each (10 records).
        for i in range(5):
            records.append(_rec(
                f"CU{i}", f"Curcgen{i} sp.",
                ["Coleoptera", "Curculionidae"]))
        for i in range(5):
            records.append(_rec(
                f"CH{i}", f"Chrygen{i} sp.",
                ["Coleoptera", "Chrysomelidae"]))
        # Diptera: 1 family, 5 genera (5 records).
        for i in range(5):
            records.append(_rec(
                f"CU2{i}", f"Culgen{i} sp.",
                ["Diptera", "Culicidae"]))
        return records

    def test_descends_through_family_to_genus(self):
        records = self._build_pool()
        items = [(r, bp._rel_lineage(r, None)) for r in records]

        chosen = bp.select(items, 6, set(), set())

        self.assertEqual(len(chosen), 6)
        orders = {r["lineage"][0] for r in chosen}
        self.assertEqual(orders, {"Coleoptera", "Diptera"})
        coleoptera_families = {
            r["lineage"][1] for r in chosen if r["lineage"][0] == "Coleoptera"
        }
        self.assertEqual(
            coleoptera_families, {"Curculionidae", "Chrysomelidae"})
        # No genus repeated.
        genera = [bp.genus_of(r["organism"]) for r in chosen]
        self.assertEqual(len(genera), len(set(genera)))


class TestOnePerGenusAtLeaf(unittest.TestCase):
    """§6.3 — at most one record per genus; longest translation wins,
    accession breaks ties."""

    def test_longest_translation_wins(self):
        records = [
            _rec("A1", "Formica rufa", ["Insecta"], translation="M" * 50),
            _rec("A2", "Formica sanguinea", ["Insecta"],
                 translation="M" * 200),
            _rec("B1", "Camponotus japonicus", ["Insecta"],
                 translation="M" * 80),
        ]
        items = [(r, bp._rel_lineage(r, None)) for r in records]

        chosen = bp.select(items, 2, set(), set())

        self.assertEqual(len(chosen), 2)
        accessions = {r["accession"] for r in chosen}
        self.assertIn("A2", accessions)   # Formica: longest wins
        self.assertNotIn("A1", accessions)
        self.assertIn("B1", accessions)

    def test_accession_breaks_length_tie(self):
        records = [
            _rec("Z9", "Formica rufa", ["Insecta"], translation="M" * 100),
            _rec("A1", "Formica sanguinea", ["Insecta"],
                 translation="M" * 100),
        ]
        items = [(r, bp._rel_lineage(r, None)) for r in records]

        chosen = bp.select(items, 1, set(), set())

        self.assertEqual([r["accession"] for r in chosen], ["A1"])

    def test_duplicate_accession_equal_length_breaks_tie_on_sequence(self):
        # Real RefSeq data has duplicate (accession, gene) CDS records
        # (e.g. isoform reannotations) with equal-length but
        # non-identical translations — accession alone cannot
        # disambiguate, so the tie-break must also be deterministic
        # under input-order shuffling (caught on the first live build,
        # task 48: `animal_mt/ATP6` picked a different one of two
        # equal-length NC_x.1 ATP6 translations depending on jsonl
        # line order).
        dup_a = _rec("D1", "Formica rufa", ["Insecta"], translation="AAAA")
        dup_b = _rec("D1", "Formica rufa", ["Insecta"], translation="AAAB")

        forward = bp.select(
            [(r, bp._rel_lineage(r, None)) for r in [dup_a, dup_b]],
            1, set(), set())
        backward = bp.select(
            [(r, bp._rel_lineage(r, None)) for r in [dup_b, dup_a]],
            1, set(), set())

        self.assertEqual(
            [r["translation"] for r in forward],
            [r["translation"] for r in backward],
        )

    def test_select_representatives_dedupes_duplicate_accession(self):
        # select_representatives() used to build its working pool as
        # `{rec["accession"]: rec for rec in candidates}` — a plain
        # dict comprehension silently keeps whichever duplicate
        # appears *last* in `candidates`, before any tie-break logic
        # ever runs, so the chosen translation depended on input
        # order. Reproduces the real failure: shuffling the input
        # order of a duplicate-accession pair (amid unrelated decoy
        # records so genus-level competition is realistic) must not
        # change which translation is kept.
        dup_a = _rec("D1", "Formica rufa", ["Insecta"], translation="AAAA")
        dup_b = _rec("D1", "Formica rufa", ["Insecta"], translation="AAAB")
        decoys = [
            _rec(f"G{i}", f"Genus{i} sp.", ["Insecta"]) for i in range(20)
        ]

        forward, _ = bp.select_representatives(
            [dup_a, dup_b] + decoys, [(None, 21)], min_per_locus=1)
        backward, _ = bp.select_representatives(
            [dup_b, dup_a] + decoys, [(None, 21)], min_per_locus=1)

        forward_d1 = [r for r in forward if r["accession"] == "D1"]
        backward_d1 = [r for r in backward if r["accession"] == "D1"]
        self.assertEqual(len(forward_d1), 1)
        self.assertEqual(
            forward_d1[0]["translation"], backward_d1[0]["translation"])


class TestUnderFill(unittest.TestCase):
    """§6.4 — budget 200 against 30 available records returns all 30,
    `selected < budget`, no exception, no padding, no duplicates."""

    def test_returns_all_available_no_padding(self):
        records = [
            _rec(f"G{i}", f"Genus{i} sp.", ["Insecta"]) for i in range(30)
        ]
        reps, stats = bp.select_representatives(
            records, [(None, 200)], min_per_locus=5)

        self.assertEqual(len(reps), 30)
        self.assertEqual(len({r["accession"] for r in reps}), 30)
        self.assertEqual(stats[0]["selected"], 30)
        self.assertLess(stats[0]["selected"], stats[0]["budget"])


class TestUnspentRedistribution(unittest.TestCase):
    """§6.5 — a sub-group that under-fills releases its unspent budget
    to a sibling that still has candidates; the loop terminates."""

    def test_shortfall_reallocated_and_terminates(self):
        # Group "A": only 2 distinct genera available.
        group_a = [
            _rec("A1", "Alpha one", ["A"]),
            _rec("A2", "Beta two", ["A"]),
        ]
        # Group "B": 50 distinct genera available.
        group_b = [
            _rec(f"B{i}", f"Genus{i:03d} sp.", ["B"]) for i in range(50)
        ]
        items = [
            (r, bp._rel_lineage(r, None)) for r in group_a + group_b
        ]

        chosen = bp.select(items, 10, set(), set())

        self.assertEqual(len(chosen), 10)
        from_a = [r for r in chosen if r["lineage"] == ["A"]]
        from_b = [r for r in chosen if r["lineage"] == ["B"]]
        self.assertEqual(len(from_a), 2)   # A's full (shortfall) capacity
        self.assertEqual(len(from_b), 8)   # picks up A's unspent budget


class TestStratumDisjointness(unittest.TestCase):
    """§6.6 — a record matching an earlier stratum's lineage term is
    never eligible for a later one; `Metazoa` after `Arthropoda` yields
    only non-arthropods."""

    def test_metazoa_after_arthropoda_excludes_arthropods(self):
        records = [
            _rec(f"AR{i}", f"Argen{i} sp.",
                 ["Metazoa", "Arthropoda", "Insecta"])
            for i in range(5)
        ] + [
            _rec(f"MO{i}", f"Molgen{i} sp.",
                 ["Metazoa", "Mollusca", "Gastropoda"])
            for i in range(5)
        ]

        chosen, stats = bp.select_representatives(
            records,
            [("Arthropoda", 3), ("Metazoa", 10)],
            min_per_locus=1,
        )

        arthropoda_stat, metazoa_stat = stats
        self.assertEqual(arthropoda_stat["selected"], 3)
        # The Metazoa stratum can only ever see the 5 mollusc records —
        # all 5 Arthropoda-matching records were removed as ineligible,
        # not just the 3 actually chosen.
        self.assertEqual(metazoa_stat["candidates"], 5)
        self.assertEqual(metazoa_stat["selected"], 5)
        for rec in chosen:
            if rec["accession"].startswith("MO"):
                self.assertNotIn("Arthropoda", rec["lineage"])


class TestDeterminism(unittest.TestCase):
    """§6.7 — same input, shuffled record order and dict construction,
    yields byte-identical (here: identical accession-ordered) output."""

    def test_shuffled_order_same_result(self):
        random.seed(12345)
        records = [
            _rec(
                f"AR{i:03d}", f"Argen{i} sp.",
                ["Metazoa", "Arthropoda", f"Order{i % 4}"])
            for i in range(40)
        ] + [
            _rec(
                f"MO{i:03d}", f"Molgen{i} sp.",
                ["Metazoa", "Mollusca"])
            for i in range(10)
        ]
        strata = [("Arthropoda", 15), ("Metazoa", 5)]

        baseline, baseline_stats = bp.select_representatives(
            list(records), strata, min_per_locus=1)

        for _ in range(5):
            shuffled = list(records)
            random.shuffle(shuffled)
            # Rebuild as dicts-of-dicts and back, to vary construction
            # order the way real jsonl-loading code would.
            by_acc = {r["accession"]: dict(r) for r in shuffled}
            reshuffled = [by_acc[k] for k in by_acc]
            result, result_stats = bp.select_representatives(
                reshuffled, strata, min_per_locus=1)

            self.assertEqual(
                [r["accession"] for r in result],
                [r["accession"] for r in baseline],
            )
            self.assertEqual(result_stats, baseline_stats)


class TestRaggedLineages(unittest.TestCase):
    """§6.8 — records with short, empty or missing lineage land in the
    unresolved bucket, are drawn last, and never crash the traversal."""

    def test_ragged_lineages_do_not_crash_and_are_selectable(self):
        well_formed = [
            _rec(f"WF{i}", f"Wellgen{i} sp.", ["Insecta", "Coleoptera"])
            for i in range(3)
        ]
        short = [_rec("SH1", "Shortgen sp.", ["Insecta"])]
        empty = [_rec("EM1", "Emptygen sp.", [])]
        missing = {
            "accession": "MI1", "organism": "Missinggen sp.",
            "gene": "COX1", "translation": "M" * 100,
        }
        records = well_formed + short + empty + [missing]

        reps, stats = bp.select_representatives(
            records, [(None, 100)], min_per_locus=1)

        self.assertEqual(len(reps), len(records))
        accessions = {r["accession"] for r in reps}
        self.assertEqual(
            accessions, {r["accession"] for r in records})

    def test_ragged_records_drawn_last_when_budget_scarce(self):
        # One well-formed record competing for genus diversity against
        # a ragged one — with budget 1, the deep-lineage record's own
        # class-level group should not starve because of the ragged
        # unresolved bucket sharing the split.
        well_formed = _rec("WF1", "Wellgen sp.", ["Insecta", "Coleoptera"])
        ragged = _rec("RG1", "Raggen sp.", [])
        items = [
            (r, bp._rel_lineage(r, None)) for r in [well_formed, ragged]
        ]
        chosen = bp.select(items, 1, set(), set())
        self.assertEqual(len(chosen), 1)


class TestMinPerLocusInteraction(unittest.TestCase):
    """§6.9 — a gene with 4 total candidates returns [] regardless of
    stratum configuration."""

    def test_below_floor_returns_empty_regardless_of_strata(self):
        records = [
            _rec(f"G{i}", f"Genus{i} sp.", ["Metazoa", "Arthropoda"])
            for i in range(4)
        ]
        reps, stats = bp.select_representatives(
            records,
            [("Arthropoda", 1000), ("Metazoa", 100)],
            min_per_locus=5,
        )
        self.assertEqual(reps, [])
        self.assertEqual(stats, [])


if __name__ == "__main__":
    unittest.main()
