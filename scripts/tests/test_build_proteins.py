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

Plus the task 64 label check (hits files hand-written; no BLAST).
"""

import importlib.util
import json
import random
import tempfile
import unittest
from argparse import Namespace
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


def _ids(*pairs):
    """[(accession, gene)] -> (gene_of, id_accession) keyed by seq id."""
    gene_of = {f"{a}_{g}": g for a, g in pairs}
    id_acc = {f"{a}_{g}": a for a, g in pairs}
    return gene_of, id_acc


def _hits(*rows):
    return [f"{q}\t{s}\t{b}\n" for q, s, b in rows]


class TestFlagMislabelled(unittest.TestCase):
    """Task 64 §3.1 — flag rule over hand-written hit tables."""

    def setUp(self):
        self.gene_of, self.id_acc = _ids(
            ("A1", "ycf3"), ("A2", "ycf3"), ("A3", "ycf3"),
            ("B1", "clpP"), ("B2", "clpP"))

    def _flag(self, rows):
        return bp.flag_mislabelled(
            self.gene_of, self.id_acc, _hits(*rows))

    def test_best_hit_other_gene_is_flagged_with_scores(self):
        flagged, no_hits = self._flag([
            ("A3_ycf3", "A1_ycf3", 41.2),
            ("A3_ycf3", "B1_clpP", 389.0),
            ("A3_ycf3", "B2_clpP", 300.0),
        ])
        self.assertEqual(flagged, [{
            "accession": "A3", "labelled_gene": "ycf3",
            "best_other_gene": "clpP",
            "own_bitscore": 41.2, "other_bitscore": 389.0,
        }])
        # A1, A2, B1, B2 have no hits at all.
        self.assertEqual(no_hits, 4)

    def test_equal_scores_not_flagged(self):
        flagged, _ = self._flag([
            ("A3_ycf3", "A1_ycf3", 100.0),
            ("A3_ycf3", "B1_clpP", 100.0),
        ])
        self.assertEqual(flagged, [])

    def test_no_own_hit_but_other_hit_is_flagged(self):
        flagged, _ = self._flag([("A3_ycf3", "B1_clpP", 50.0)])
        self.assertEqual(len(flagged), 1)
        self.assertIsNone(flagged[0]["own_bitscore"])

    def test_no_hits_counted_not_flagged(self):
        flagged, no_hits = self._flag([])
        self.assertEqual(flagged, [])
        self.assertEqual(no_hits, 5)

    def test_self_hits_ignored(self):
        flagged, no_hits = self._flag([("A1_ycf3", "A1_ycf3", 999.0)])
        self.assertEqual(flagged, [])
        self.assertEqual(no_hits, 5)

    def test_multiple_hsps_take_max_and_input_order_irrelevant(self):
        rows = [
            ("A3_ycf3", "A1_ycf3", 10.0),
            ("A3_ycf3", "A1_ycf3", 90.0),
            ("A3_ycf3", "B1_clpP", 80.0),
        ]
        forward, _ = self._flag(rows)
        backward, _ = self._flag(list(reversed(rows)))
        self.assertEqual(forward, [])
        self.assertEqual(forward, backward)


class TestExclusionReselection(unittest.TestCase):
    """Task 64 §3.1 — re-selection after exclusion."""

    def _pool(self):
        return [
            _rec(f"G{i}", f"Genus{i} sp.", ["Insecta"],
                 translation="M" * (100 - i))
            for i in range(6)
        ]

    def test_budget_refilled_with_next_candidate(self):
        pool = self._pool()
        strata = [(None, 3)]
        before, _ = bp.select_representatives(pool, strata, 3)
        self.assertEqual([r["accession"] for r in before],
                         ["G0", "G1", "G2"])
        trimmed = [r for r in pool if r["accession"] != "G1"]
        after, _ = bp.select_representatives(trimmed, strata, 3)
        self.assertEqual([r["accession"] for r in after],
                         ["G0", "G2", "G3"])

    def test_reselection_independent_of_input_order(self):
        trimmed = [r for r in self._pool() if r["accession"] != "G1"]
        expected, _ = bp.select_representatives(trimmed, [(None, 3)], 3)
        shuffled = list(trimmed)
        random.Random(7).shuffle(shuffled)
        got, _ = bp.select_representatives(shuffled, [(None, 3)], 3)
        self.assertEqual(expected, got)


class TestPipeline(unittest.TestCase):
    """Task 64 §3.1 — full / label_check round trip on a tiny panel."""

    GENE_SETS = {"sets": {"plant_pt": {
        "protein_coding": ["ycf3", "clpP", "rbcL"]}}}

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = bp.Path(self._tmp.name)
        (self.tmp / "gs.json").write_text(json.dumps(self.GENE_SETS))
        recs = []
        for gene, n in (("ycf3", 4), ("clpP", 4), ("rbcL", 3)):
            for i in range(n):
                recs.append(_rec(
                    f"{gene}{i}", f"{gene}gen{i} sp.", ["Viridiplantae"],
                    gene=gene, translation="M" * (50 + i)))
        # Mislabelled: really clpP, labelled ycf3.
        recs.append(_rec("BAD", "Badgen sp.", ["Viridiplantae"],
                         gene="ycf3", translation="M" * 99))
        (self.tmp / "cds.jsonl").write_text(
            "\n".join(json.dumps(r) for r in recs) + "\n")

    def _args(self, **kw):
        base = dict(
            gene_sets=str(self.tmp / "gs.json"),
            cds_jsonl=[f"plant_pt={self.tmp / 'cds.jsonl'}"],
            out=str(self.tmp / "out"), refseq_release="237",
            min_per_locus=3, max_per_locus=10, stratum=[],
            exclude="", emit_selection=str(self.tmp / "sel"))
        base.update(kw)
        return Namespace(**base)

    def _lc(self, pass_number, max_passes=3):
        return Namespace(
            selection=str(self.tmp / "sel"),
            exclusions=str(self.tmp / "ex.json"),
            pass_number=pass_number, max_passes=max_passes)

    def _write_hits(self, rows):
        (self.tmp / "sel" / "plant_pt.hits.tsv").write_text(
            "".join(_hits(*rows)))

    def test_emit_flag_exclude_final(self):
        bp.main_full(self._args())
        sel = json.loads((self.tmp / "sel" / "selection.json").read_text())
        self.assertIn(
            "BAD", [r["accession"] for r in sel["plant_pt"]["ycf3"]])
        self.assertTrue((self.tmp / "sel" / "plant_pt.faa").is_file())

        self._write_hits([
            ("BAD_ycf3", "clpP0_clpP", 300.0),
            ("BAD_ycf3", "ycf30_ycf3", 40.0),
        ])
        self.assertEqual(bp.main_label_check(self._lc(1)), bp.EXIT_FLAGGED)

        # Pass 2 is clean once BAD has been replaced.
        bp.main_full(self._args(exclude=str(self.tmp / "ex.json")))
        self._write_hits([])
        self.assertEqual(bp.main_label_check(self._lc(2)), 0)

        bp.main_full(self._args(
            exclude=str(self.tmp / "ex.json"), emit_selection=""))
        prov = json.loads((self.tmp / "out" / "provenance.json").read_text())
        lc = prov["label_check"]["plant_pt"]
        self.assertEqual(lc["iterations"], 2)
        self.assertEqual(lc["method"], bp.LABEL_CHECK_METHOD)
        self.assertEqual(
            [e["accession"] for e in lc["excluded"]], ["BAD"])
        self.assertEqual(lc["excluded"][0]["best_other_gene"], "clpP")
        self.assertNotIn(
            "BAD", prov["genes"]["plant_pt"]["ycf3"]["accessions"])
        self.assertEqual(
            prov["genes"]["plant_pt"]["ycf3"]["candidates"], 4)

    def test_gene_pushed_below_minimum_is_dropped_and_recorded(self):
        (self.tmp / "ex.json").write_text(json.dumps({"plant_pt": {
            "iterations": 1, "no_hits": 0, "excluded": [{
                "accession": "rbcL0", "labelled_gene": "rbcL",
                "best_other_gene": "clpP",
                "own_bitscore": 1.0, "other_bitscore": 2.0}]}}))
        bp.main_full(self._args(
            exclude=str(self.tmp / "ex.json"), emit_selection=""))
        prov = json.loads((self.tmp / "out" / "provenance.json").read_text())
        self.assertNotIn("rbcL", prov["genes"]["plant_pt"])
        self.assertFalse((self.tmp / "out" / "plant_pt" / "rbcL.faa").exists())
        self.assertEqual(
            prov["label_check"]["plant_pt"]["dropped_genes"],
            [{"gene": "rbcL", "candidates_before_exclusion": 3,
              "candidates_after_exclusion": 2}])

    def test_fails_at_iteration_cap(self):
        bp.main_full(self._args())
        self._write_hits([("BAD_ycf3", "clpP0_clpP", 300.0)])
        self.assertEqual(bp.main_label_check(self._lc(3, 3)), 1)
        self.assertEqual(bp.main_label_check(self._lc(2, 3)), bp.EXIT_FLAGGED)


if __name__ == "__main__":
    unittest.main()
