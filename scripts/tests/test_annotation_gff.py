"""Unit tests for bin/annotation_gff.py — task 40 §2 (shared parser).

Cases:
  1. Miniprot-provenance mRNA (Target=, no Name=) — gene from Target.
  2. Rescued (mitos)-provenance mRNA (Name=) — gene from Name.
  3. Multi-exon feature — CDS intervals accumulate onto its parent.
  4. mRNA with no CDS children — dropped from the result.
  5. CDS row with unknown Parent — ignored, no KeyError.
  6. mRNA row missing ID= — skipped.
  7. Non-9-column / comment / blank lines — skipped, not fatal.
  8. Non-numeric start/end on mRNA or CDS — skipped, not fatal.
  9. parse_attrs: entries without '=' are ignored.
"""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[2] / "bin" / "annotation_gff.py"
)
_spec = importlib.util.spec_from_file_location(
    "annotation_gff", MODULE_PATH)
agff = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(agff)


class ParseAnnotationGffTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _write(self, lines):
        path = self.tmp / "merged.gff"
        path.write_text("\n".join(lines) + "\n")
        return path

    def test_miniprot_gene_from_target(self):
        path = self._write([
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\t"
            "ID=MP1;Target=panel_ATP6 1 30",
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP1",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]["gene"], "ATP6")
        self.assertEqual(features[0]["source"], "miniprot")
        self.assertEqual(features[0]["cds"], [(1, 90)])

    def test_rescued_gene_from_name(self):
        path = self._write([
            "ctg1\tmitos\tmRNA\t1\t50\t.\t+\t.\t"
            "ID=mitos_rescue_ATP8_1;Name=ATP8;OriginalName=atp8_0",
            "ctg1\tmitos\tCDS\t1\t50\t.\t+\t0\t"
            "ID=mitos_rescue_ATP8_1_cds1;Parent=mitos_rescue_ATP8_1",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]["gene"], "ATP8")
        self.assertEqual(features[0]["source"], "mitos")

    def test_multi_exon_accumulates_cds(self):
        path = self._write([
            "ctg1\tmitos\tmRNA\t1\t100\t.\t+\t.\tID=r1;Name=NAD1",
            "ctg1\tmitos\tCDS\t1\t40\t.\t+\t0\tParent=r1",
            "ctg1\tmitos\tCDS\t60\t100\t.\t+\t1\tParent=r1",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(features[0]["cds"], [(1, 40), (60, 100)])

    def test_mrna_without_cds_dropped(self):
        path = self._write([
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\t"
            "ID=MP1;Target=panel_ATP6 1 30",
        ])
        self.assertEqual(agff.parse_annotation_gff(path), [])

    def test_cds_with_unknown_parent_ignored(self):
        path = self._write([
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=ghost",
        ])
        self.assertEqual(agff.parse_annotation_gff(path), [])

    def test_mrna_missing_id_skipped(self):
        path = self._write([
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\tTarget=panel_ATP6 1 30",
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP1",
        ])
        self.assertEqual(agff.parse_annotation_gff(path), [])

    def test_malformed_lines_skipped(self):
        path = self._write([
            "##gff-version 3",
            "",
            "ctg1\tminiprot\tmRNA\ttoo\tfew\tcolumns",
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\t"
            "ID=MP1;Target=panel_ATP6 1 30",
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP1",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(len(features), 1)

    def test_non_numeric_coordinates_skipped(self):
        path = self._write([
            "ctg1\tminiprot\tmRNA\tNA\tNA\t.\t+\t.\t"
            "ID=MP1;Target=panel_ATP6 1 30",
            "ctg1\tminiprot\tmRNA\t1\t90\t.\t+\t.\t"
            "ID=MP2;Target=panel_ATP6 1 30",
            "ctg1\tminiprot\tCDS\tNA\tNA\t.\t+\t0\tParent=MP2",
            "ctg1\tminiprot\tCDS\t1\t90\t.\t+\t0\tParent=MP2",
        ])
        features = agff.parse_annotation_gff(path)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]["id"], "MP2")
        self.assertEqual(features[0]["cds"], [(1, 90)])

    def test_parse_attrs_ignores_entries_without_equals(self):
        attrs = agff.parse_attrs("ID=x;standalone;Name=y")
        self.assertEqual(attrs, {"ID": "x", "Name": "y"})

    def test_gene_symbol_no_target_no_name(self):
        self.assertEqual(agff.gene_symbol({}), "")


def rec(gene, seqid, strand, exons, identity, mrna_id=""):
    return {
        "gene": gene, "seqid": seqid, "strand": strand, "cds": exons,
        "identity": identity, "mrna_id": mrna_id,
    }


def cluster(records):
    return agff.cluster_loci(records, lambda r: r["gene"])


def ids(winners):
    return sorted(w["mrna_id"] for w in winners)


class ClusterLociTests(unittest.TestCase):

    def test_rpl2_bridging_hit_does_not_swallow_other_copy(self):
        # Two IR copies on opposite strands plus a hit whose last exon
        # sits 44 kb upstream, so its outer span covers both copies.
        plus = [(109894, 110283), (110937, 111380)]
        minus = [(129890, 130333), (130987, 131376)]
        bridge = [(130987, 131376), (129901, 130333), (85769, 85881)]
        for bridge_identity, minus_winner in ((0.9455, "minus"),
                                              (0.99, "bridge")):
            winners = cluster([
                rec("rpl2", "p", "+", plus, 0.958, "plus"),
                rec("rpl2", "p", "-", minus, 0.95, "minus"),
                rec("rpl2", "p", "-", bridge, bridge_identity, "bridge"),
            ])["rpl2"]
            self.assertEqual(ids(winners), sorted(["plus", minus_winner]))

    def test_opposite_strands_are_separate_loci(self):
        winners = cluster([
            rec("g", "c", "+", [(1, 100)], 0.9, "a"),
            rec("g", "c", "-", [(1, 100)], 0.8, "b"),
        ])["g"]
        self.assertEqual(ids(winners), ["a", "b"])

    def test_hit_inside_an_intron_is_a_separate_locus(self):
        winners = cluster([
            rec("g", "c", "+", [(1, 100), (500, 600)], 0.9, "outer"),
            rec("g", "c", "+", [(200, 300)], 0.8, "inner"),
        ])["g"]
        self.assertEqual(ids(winners), ["inner", "outer"])

    def test_different_contigs_are_separate_loci(self):
        winners = cluster([
            rec("g", "c1", "+", [(1, 100)], 0.9, "a"),
            rec("g", "c2", "+", [(1, 100)], 0.8, "b"),
        ])["g"]
        self.assertEqual(ids(winners), ["a", "b"])

    def test_exon_overlap_merges_and_highest_identity_wins(self):
        winners = cluster([
            rec("g", "c", "+", [(1, 100), (500, 600)], 0.7, "a"),
            rec("g", "c", "+", [(5, 100), (505, 600)], 0.9, "b"),
            rec("g", "c", "+", [(1, 98), (500, 598)], 0.8, "c"),
        ])["g"]
        self.assertEqual(ids(winners), ["b"])

    def test_merge_is_transitive(self):
        winners = cluster([
            rec("g", "c", "+", [(1, 10)], 0.5, "a"),
            rec("g", "c", "+", [(8, 20)], 0.6, "b"),
            rec("g", "c", "+", [(15, 30)], 0.7, "c"),
        ])["g"]
        self.assertEqual(ids(winners), ["c"])

    def test_adjacent_but_disjoint_exons_stay_separate(self):
        winners = cluster([
            rec("g", "c", "+", [(1, 10)], 0.5, "a"),
            rec("g", "c", "+", [(11, 20)], 0.6, "b"),
        ])["g"]
        self.assertEqual(ids(winners), ["a", "b"])

    def test_tie_is_independent_of_input_order(self):
        hits = [
            rec("g", "c", "+", [(10, 100)], 0.9, "late"),
            rec("g", "c", "+", [(1, 100)], 0.9, "zzz"),
            rec("g", "c", "+", [(1, 100)], 0.9, "aaa"),
        ]
        forward = cluster([dict(h) for h in hits])["g"]
        backward = cluster([dict(h) for h in reversed(hits)])["g"]
        self.assertEqual(ids(forward), ["aaa"])
        self.assertEqual(ids(backward), ["aaa"])

    def test_same_gene_in_two_places_is_two_features(self):
        winners = cluster([
            rec("ndhB", "c", "+", [(100, 200)], 0.9, "a"),
            rec("ndhB", "c", "+", [(9000, 9100)], 0.9, "b"),
        ])["ndhB"]
        self.assertEqual(len(winners), 2)

    def test_sets_start_and_end(self):
        winner = cluster([
            rec("g", "c", "+", [(50, 60), (10, 20)], 0.9, "a"),
        ])["g"][0]
        self.assertEqual((winner["start"], winner["end"]), (10, 60))

    def test_missing_mrna_id_is_tolerated_in_tie(self):
        a = {"gene": "g", "seqid": "c", "strand": "+", "cds": [(1, 9)],
             "identity": 0.5}
        b = dict(a)
        self.assertEqual(len(cluster([a, b])["g"]), 1)

    def test_exons_overlap_predicate(self):
        self.assertTrue(agff.exons_overlap([(1, 10)], [(10, 20)]))
        self.assertFalse(agff.exons_overlap([(1, 10)], [(11, 20)]))
        self.assertFalse(
            agff.exons_overlap([(1, 10), (50, 60)], [(20, 30)]))

    def test_sweep_agrees_with_pairwise_predicate(self):
        import random
        rng = random.Random(7)
        for _ in range(50):
            recs = []
            for n in range(8):
                starts = sorted(rng.sample(range(0, 200, 5), 2))
                first_end = min(starts[0] + rng.randint(0, 12),
                                starts[1] - 1)
                exons = [(starts[0], first_end),
                         (starts[1], starts[1] + rng.randint(0, 12))]
                recs.append(rec("g", "c", "+", exons, 0.5, str(n)))
            comps = agff._components(recs)
            where = {r["mrna_id"]: i for i, c in enumerate(comps)
                     for r in c}
            for a in recs:
                for b in recs:
                    if agff.exons_overlap(a["cds"], b["cds"]):
                        self.assertEqual(
                            where[a["mrna_id"]], where[b["mrna_id"]])
            for comp in comps:
                seen = {comp[0]["mrna_id"]}
                grew = True
                while grew:
                    grew = False
                    for a in comp:
                        for b in comp:
                            if (a["mrna_id"] in seen
                                    and b["mrna_id"] not in seen
                                    and agff.exons_overlap(
                                        a["cds"], b["cds"])):
                                seen.add(b["mrna_id"])
                                grew = True
                self.assertEqual(len(seen), len(comp))


RPL2_GFF = "\n".join([
    "##gff-version 3",
    "p\tminiprot\tmRNA\t109894\t111380\t100\t+\t.\t"
    "ID=MP1;Identity=0.9580;Target=A_rpl2 1 50",
    "p\tminiprot\tCDS\t109894\t110283\t100\t+\t0\tParent=MP1",
    "p\tminiprot\tCDS\t110937\t111380\t100\t+\t0\tParent=MP1",
    "p\tminiprot\tmRNA\t85769\t131376\t100\t-\t.\t"
    "ID=MP2;Identity=0.9990;Target=B_rpl2 1 50",
    "p\tminiprot\tCDS\t130987\t131376\t100\t-\t0\tParent=MP2",
    "p\tminiprot\tCDS\t129901\t130333\t100\t-\t0\tParent=MP2",
    "p\tminiprot\tCDS\t85769\t85881\t100\t-\t0\tParent=MP2",
    "p\tminiprot\tmRNA\t129890\t131376\t100\t-\t.\t"
    "ID=MP3;Identity=0.9455;Target=C_rpl2 1 50",
    "p\tminiprot\tCDS\t130987\t131376\t100\t-\t0\tParent=MP3",
    "p\tminiprot\tCDS\t129890\t130333\t100\t-\t0\tParent=MP3",
]) + "\n"


class BarcodeAnnotationAgreementTests(unittest.TestCase):
    """C5 and C8 pick the same winners from the same cds.gff."""

    def test_same_winners_from_both_callers(self):
        sys.path.insert(0, str(MODULE_PATH.parent))
        import annotate_summary as asum
        import validate_barcodes as vb

        gff = Path(tempfile.mkdtemp()) / "cds.gff"
        gff.write_text(RPL2_GFF)
        c5 = vb.cluster_by_locus(vb.parse_gff(gff))["RPL2"]
        c8 = asum.cluster_cds_by_gene(asum.parse_cds_gff(gff))
        c8 = [w for ws in c8.values() for w in ws]

        def key(w):
            return (w["mrna_id"], w["start"], w["end"])

        self.assertEqual(sorted(map(key, c5)), sorted(map(key, c8)))
        self.assertEqual(sorted(w["mrna_id"] for w in c5), ["MP1", "MP2"])


if __name__ == "__main__":
    unittest.main()
