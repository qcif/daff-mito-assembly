"""Unit tests for bin/rotate_target.py — spec §2.2 C13, §3.8, tasks 75, 76.

Fixtures are small synthetic FASTA / miniprot-shaped GFF3 fragments
written by the tests, not reference data (CONSTITUTION rule 19). The
tandem-duplicate tests also use two real contigs kept under
fixtures/tandem_duplicate/.
"""

import contextlib
import importlib.util
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

MODULE_PATH = (
    Path(__file__).resolve().parents[2] / "bin" / "rotate_target.py"
)
_spec = importlib.util.spec_from_file_location("rotate_target", MODULE_PATH)
rt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rt)

ANCHORS = ["COX1", "COX2", "COX3", "CYTB"]
CONTIG = "contig_10"
HEADER = ">contig_10 length=1000 circ=Y"
LENGTH = 1000
FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "tandem_duplicate"
TANDEM = {
    "kmer": 31, "min_support": 0.5, "tolerance_bp": 50,
    "min_period_bp": 1000,
}
GENOME_LEN = 15000


def make_seq(length=LENGTH, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def gff_line(
    gene, start, end, strand="+", score=500, query_start=1, contig=CONTIG,
):
    attrs = f"ID=MP000001;Rank=1;Target=NC_000000_{gene} {query_start} 500"
    return "\t".join([
        contig, "miniprot", "mRNA", str(start), str(end), str(score),
        strand, ".", attrs,
    ])


def circular_meta(circular=True, method="flye_circ", cid=CONTIG):
    return {"contigs": [{
        "contig_id": cid, "circular": circular, "circular_method": method,
    }]}


def hit(gene="COX1", **kwargs):
    base = {
        "contig": CONTIG, "gene": gene, "table": 5, "score": 500.0,
        "start": 100, "end": 400, "strand": "+", "query_start": 1,
    }
    base.update(kwargs)
    return base


def is_cyclic_shift(a, b):
    return len(a) == len(b) and b in a + a


class RotateAllMixin:
    def rotate(self, seq, hits, anchors=ANCHORS, meta=None, header=HEADER):
        meta = circular_meta() if meta is None else meta
        out, rotation = rt.rotate_all(
            [(header, seq)], meta, hits, anchors)
        return out[0][1], rotation["contigs"][0], rotation


class TestRotation(RotateAllMixin, unittest.TestCase):
    def test_plus_strand_anchor_starts_at_one(self):
        seq = make_seq()
        out, rec, rotation = self.rotate(seq, [hit(start=101, end=400)])
        self.assertEqual(out, seq[100:] + seq[:100])
        self.assertEqual(rec["status"], "rotated")
        self.assertEqual(rec["offset"], 100)
        self.assertEqual(rec["anchor_gene"], "COX1")
        self.assertFalse(rec["fallback"])
        self.assertFalse(rec["reverse_complemented"])
        self.assertEqual(rec["anchor_table"], 5)
        self.assertEqual(rec["length_bp"], LENGTH)
        self.assertEqual(rotation["status"], "rotated")

    def test_minus_strand_anchor_is_reverse_complemented(self):
        seq = make_seq()
        # On the minus strand the CDS starts at its end coordinate.
        out, rec, _ = self.rotate(
            seq, [hit(start=101, end=400, strand="-")])
        flipped = rt.revcomp(seq)
        self.assertEqual(out, flipped[600:] + flipped[:600])
        self.assertTrue(rec["reverse_complemented"])
        self.assertEqual(rec["raw_start"], 601)
        self.assertEqual(rec["offset"], 600)
        # The anchor's first codon is now at position 1 on the + strand.
        self.assertEqual(out[:3], rt.revcomp(seq[397:400]))

    def test_anchor_already_at_position_one(self):
        seq = make_seq()
        out, rec, _ = self.rotate(seq, [hit(start=1, end=300)])
        self.assertEqual(out, seq)
        self.assertEqual(rec["status"], "rotated")
        self.assertEqual(rec["offset"], 0)

    def test_query_start_extends_upstream(self):
        seq = make_seq()
        out, rec, _ = self.rotate(
            seq, [hit(start=106, end=400, query_start=4)])
        self.assertEqual(rec["raw_start"], 106)
        self.assertEqual(rec["query_start"], 4)
        self.assertEqual(rec["offset"], 96)
        self.assertEqual(out, seq[96:] + seq[:96])

    def test_query_start_extension_wraps_across_origin(self):
        seq = make_seq()
        out, rec, _ = self.rotate(
            seq, [hit(start=4, end=300, query_start=3)])
        self.assertEqual(rec["offset"], LENGTH - 3)
        self.assertEqual(out, seq[-3:] + seq[:-3])

    def test_minus_strand_with_query_start(self):
        seq = make_seq()
        out, rec, _ = self.rotate(
            seq, [hit(start=101, end=400, strand="-", query_start=2)])
        flipped = rt.revcomp(seq)
        self.assertEqual(rec["offset"], 597)
        self.assertEqual(out, flipped[597:] + flipped[:597])


class TestAnchorChoice(RotateAllMixin, unittest.TestCase):
    def test_fallback_when_cox1_absent(self):
        seq = make_seq()
        _, rec, _ = self.rotate(seq, [hit("COX2", start=201, end=500)])
        self.assertEqual(rec["anchor_gene"], "COX2")
        self.assertTrue(rec["fallback"])

    def test_order_beats_score(self):
        seq = make_seq()
        hits = [
            hit("COX2", start=201, end=500, score=900),
            hit("COX1", start=601, end=900, score=100),
        ]
        _, rec, _ = self.rotate(seq, hits)
        self.assertEqual(rec["anchor_gene"], "COX1")
        self.assertEqual(rec["offset"], 600)

    def test_highest_score_wins_then_start_then_table(self):
        hits = [
            hit(start=50, score=100, table=5),
            hit(start=60, score=300, table=5),
            hit(start=70, score=300, table=2),
        ]
        self.assertEqual(rt.best_hit(hits)["start"], 60)
        tied = [hit(start=80, table=5), hit(start=70, table=5)]
        self.assertEqual(rt.best_hit(tied)["start"], 70)
        same = [hit(start=70, table=5), hit(start=70, table=2)]
        self.assertEqual(rt.best_hit(same)["table"], 2)
        self.assertIsNone(rt.best_hit([]))

    def test_non_anchor_gene_ignored(self):
        _, rec, _ = self.rotate(make_seq(), [hit("ND1")])
        self.assertEqual(rec["status"], "no_anchor")


class TestPassthroughStatuses(RotateAllMixin, unittest.TestCase):
    def test_no_hits_is_no_anchor(self):
        seq = make_seq()
        out, rec, rotation = self.rotate(seq, [])
        self.assertEqual(out, seq)
        self.assertEqual(rec["status"], "no_anchor")
        self.assertEqual(rec["offset"], 0)
        self.assertIsNone(rec["anchor_gene"])
        self.assertEqual(rotation["status"], "no_anchor")

    def test_linear_contig_is_not_circular(self):
        seq = make_seq()
        for meta in (
            circular_meta(circular=False, method="none"),
            circular_meta(circular=True, method="none"),
            {"contigs": []},
        ):
            out, rec, _ = self.rotate(seq, [hit()], meta=meta)
            self.assertEqual(out, seq)
            self.assertEqual(rec["status"], "not_circular")

    def test_end_overlap_counts_as_circular(self):
        _, rec, _ = self.rotate(
            make_seq(), [hit(start=11)],
            meta=circular_meta(method="end_overlap"))
        self.assertEqual(rec["status"], "rotated")

    def test_empty_fasta_is_no_target(self):
        out, rotation = rt.rotate_all([], circular_meta(), [], ANCHORS)
        self.assertEqual(out, [])
        self.assertEqual(rotation["status"], "no_target")
        self.assertEqual(rotation["contigs"], [])

    def test_empty_anchor_list_is_not_configured(self):
        seq = make_seq()
        out, rec, rotation = self.rotate(seq, [hit()], anchors=[])
        self.assertEqual(out, seq)
        self.assertEqual(rec["status"], "not_configured")
        self.assertEqual(rotation["status"], "not_configured")

    def test_empty_sequence_record_is_no_anchor(self):
        record, seq = rt.process_contig(
            ">c", "", circular_meta(cid="c")["contigs"][0],
            [hit(contig="c")], ANCHORS)
        self.assertEqual(record["status"], "no_anchor")
        self.assertEqual(seq, "")


class TestMultiContig(RotateAllMixin, unittest.TestCase):
    def test_contigs_handled_independently(self):
        seq_a, seq_b = make_seq(1000, 1), make_seq(1200, 2)
        meta = {"contigs": [
            {"contig_id": "contig_10", "circular": True,
             "circular_method": "flye_circ"},
            {"contig_id": "contig_11", "circular": True,
             "circular_method": "flye_circ"},
        ]}
        hits = [
            hit(contig="contig_10", start=101, end=400),
            hit(contig="absent", start=5, end=50),
        ]
        out, rotation = rt.rotate_all(
            [(">contig_10 x", seq_a), (">contig_11", seq_b)],
            meta, hits, ANCHORS)
        self.assertEqual(out[0], (">contig_10 x", seq_a[100:] + seq_a[:100]))
        self.assertEqual(out[1], (">contig_11", seq_b))
        statuses = [c["status"] for c in rotation["contigs"]]
        self.assertEqual(statuses, ["rotated", "no_anchor"])
        # Primary (longest) contig decides the top-level status.
        self.assertEqual(rotation["status"], "no_anchor")

    def test_invariants(self):
        rng = random.Random(7)
        for _ in range(50):
            length = rng.randint(300, 1500)
            seq = make_seq(length, rng.randint(0, 10**6))
            start = rng.randint(1, length - 100)
            strand = rng.choice("+-")
            h = hit(
                start=start, end=start + 90, strand=strand,
                query_start=rng.randint(1, 20))
            out, rec, _ = self.rotate(seq, [h])
            base = rt.revcomp(seq) if strand == "-" else seq
            self.assertEqual(len(out), length)
            self.assertEqual(sorted(out), sorted(base))
            self.assertTrue(is_cyclic_shift(base, out))
            self.assertEqual(rec["length_bp"], length)


class TestParsing(unittest.TestCase):
    def test_parse_hits(self):
        lines = [
            "##gff-version 3",
            "",
            gff_line("COX1", 10, 300, "-", 700, 3),
            gff_line("ND1", 10, 300),
            "too\tfew\tcolumns",
            "\t".join([CONTIG, "miniprot", "CDS", "1", "2", "3", "+", ".",
                       "Parent=MP1"]),
            gff_line("COX2", 10, 300).replace("Target=", "Other="),
            gff_line("COX2", 10, 300).replace(" 1 500", ""),
            gff_line("COX2", 10, 300).replace("\t500\t", "\tNA\t"),
            gff_line("COX2", 10, 300, strand="."),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.gff"
            path.write_text("\n".join(lines) + "\n")
            hits = rt.parse_hits(path, 5, ANCHORS)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene"], "COX1")
        self.assertEqual(hits[0]["strand"], "-")
        self.assertEqual(hits[0]["query_start"], 3)
        self.assertEqual(hits[0]["table"], 5)
        self.assertEqual(hits[0]["score"], 700.0)

    def test_read_write_fasta_preserves_header(self):
        seq = make_seq(205)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "t.fa"
            path.write_text("stray\n>a desc  here\n" + seq + "\n>b\n\n")
            records = rt.read_fasta(path)
            self.assertEqual(records, [(">a desc  here", seq), (">b", "")])
            rt.write_fasta(path, records)
            self.assertEqual(rt.read_fasta(path), records)
            self.assertEqual(rt.contig_id_of(">a desc"), "a")
            self.assertEqual(rt.contig_id_of(">"), "")

    def test_read_empty_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "t.fa"
            path.write_text("")
            self.assertEqual(rt.read_fasta(path), [])

    def test_parse_gff_arg(self):
        self.assertEqual(rt.parse_gff_arg("5:a/b.gff"), (5, Path("a/b.gff")))


def second_copy(genome, seed):
    """A diverged copy whose first k-mer is intact."""
    return genome[:31] + mutate(genome[31:], 0.005, seed)


def mutate(seq, rate, seed):
    rng = random.Random(seed)
    return "".join(
        rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b
        for b in seq)


class TandemMixin:
    def run_tandem(
        self, seq, hits=(), anchors=ANCHORS, params=TANDEM, meta=None,
    ):
        meta = circular_meta() if meta is None else meta
        out, rotation = rt.rotate_all(
            [(HEADER, seq)], meta, list(hits), anchors, params)
        return out[0][1], rotation["contigs"][0]


class TestTandemSynthetic(TandemMixin, unittest.TestCase):
    def setUp(self):
        self.genome = make_seq(GENOME_LEN, seed=7)

    def test_two_copies_collapse_to_one(self):
        seq = self.genome + second_copy(self.genome, seed=1)
        out, rec = self.run_tandem(seq, anchors=[])
        tandem = rec["tandem"]
        self.assertEqual(tandem["status"], "collapsed")
        self.assertEqual(tandem["copies"], 2)
        self.assertAlmostEqual(tandem["period_bp"], GENOME_LEN, delta=5)
        self.assertEqual(tandem["cut_method"], "leading_kmer")
        self.assertEqual(out, seq[:tandem["cut_bp"]])
        self.assertEqual(rec["length_bp"], len(out))
        self.assertEqual(rec["input_length_bp"], len(seq))
        self.assertEqual(tandem["params"], TANDEM)

    def test_three_copies_keep_one(self):
        seq = (
            self.genome + second_copy(self.genome, seed=1)
            + second_copy(self.genome, seed=2))
        out, rec = self.run_tandem(seq, anchors=[])
        self.assertEqual(rec["tandem"]["status"], "collapsed")
        self.assertEqual(rec["tandem"]["copies"], 3)
        self.assertAlmostEqual(len(out), GENOME_LEN, delta=5)

    def test_segmental_duplication_is_single(self):
        genome = make_seq(30000, seed=8)
        seq = genome + genome[5000:10000]
        out, rec = self.run_tandem(seq, anchors=[])
        self.assertEqual(rec["tandem"]["status"], "single")
        self.assertEqual(out, seq)
        self.assertEqual(rec["length_bp"], len(seq))
        self.assertIsNone(rec["tandem"]["cut_bp"])

    def test_non_tiling_period_is_unresolved(self):
        seq = self.genome + self.genome + make_seq(3000, seed=9)
        out, rec = self.run_tandem(seq, anchors=[])
        self.assertEqual(rec["tandem"]["status"], "unresolved")
        self.assertEqual(out, seq)
        self.assertEqual(rec["length_bp"], len(seq))

    def test_microsatellite_is_single(self):
        seq = make_seq(5000, seed=10) + "AT" * 3000
        out, rec = self.run_tandem(seq, anchors=[])
        self.assertEqual(rec["tandem"]["status"], "single")
        self.assertEqual(out, seq)

    def test_no_repeats_has_no_period(self):
        out, rec = self.run_tandem(make_seq(5000, seed=11), anchors=[])
        self.assertEqual(rec["tandem"]["status"], "single")
        self.assertIsNone(rec["tandem"]["period_bp"])
        self.assertIsNone(rec["tandem"]["copies"])

    def test_support_equal_to_threshold_is_repetitive(self):
        status = rt.tandem_decision(2000, 1000, 0.5, 2, TANDEM)
        self.assertEqual(status, "collapsed")
        status = rt.tandem_decision(2000, 1000, 0.4999, 2, TANDEM)
        self.assertEqual(status, "single")
        status = rt.tandem_decision(2000, 999, 0.9, 2, TANDEM)
        self.assertEqual(status, "single")

    def test_one_copy_period_is_unresolved(self):
        status = rt.tandem_decision(1500, 1400, 0.9, 1, TANDEM)
        self.assertEqual(status, "unresolved")

    def test_mutated_leading_kmer_falls_back_to_period(self):
        copy2 = (
            mutate(self.genome[:31], 1.0, seed=3) + self.genome[31:])
        seq = self.genome + copy2
        out, rec = self.run_tandem(seq, anchors=[])
        self.assertEqual(rec["tandem"]["status"], "collapsed")
        self.assertEqual(rec["tandem"]["cut_method"], "period")
        self.assertEqual(len(out), rec["tandem"]["period_bp"])

    def test_leading_kmer_far_from_period_falls_back(self):
        copy2 = (
            mutate(self.genome[:31], 1.0, seed=3) + self.genome[31:])
        seq = self.genome + copy2 + self.genome
        out, rec = self.run_tandem(seq, anchors=[])
        self.assertEqual(rec["tandem"]["status"], "collapsed")
        self.assertEqual(rec["tandem"]["cut_method"], "period")
        self.assertAlmostEqual(len(out), GENOME_LEN, delta=5)

    def test_not_tested_cases(self):
        seq = self.genome + self.genome
        out, rec = self.run_tandem(
            seq, anchors=[], meta=circular_meta(circular=False))
        self.assertEqual(rec["tandem"]["status"], "not_tested")
        self.assertEqual(out, seq)
        out, rec = self.run_tandem(seq, anchors=[], params=None)
        self.assertEqual(rec["tandem"]["status"], "not_tested")
        self.assertIsNone(rec["tandem"]["params"])
        self.assertEqual(out, seq)
        out, rec = self.run_tandem("", anchors=[])
        self.assertEqual(rec["tandem"]["status"], "not_tested")
        self.assertEqual(out, "")

    def test_no_anchor_contig_collapses_in_flye_frame(self):
        seq = self.genome + second_copy(self.genome, seed=1)
        out, rec = self.run_tandem(seq)
        self.assertEqual(rec["status"], "no_anchor")
        self.assertEqual(rec["tandem"]["status"], "collapsed")
        self.assertEqual(out, seq[:rec["tandem"]["cut_bp"]])

    def test_anchor_in_second_copy_starts_kept_copy(self):
        seq = self.genome + second_copy(self.genome, seed=1)
        start = GENOME_LEN + 2001
        out, rec = self.run_tandem(
            seq, [hit(start=start, end=start + 299)])
        self.assertEqual(rec["status"], "rotated")
        self.assertEqual(rec["tandem"]["status"], "collapsed")
        self.assertEqual(rec["tandem"]["cut_method"], "leading_kmer")
        self.assertTrue(out.startswith(seq[start - 1:start + 30]))
        self.assertAlmostEqual(len(out), GENOME_LEN, delta=5)
        self.assertEqual(rec["length_bp"], len(out))
        self.assertEqual(rec["input_length_bp"], len(seq))


class TestTandemFixtures(TandemMixin, unittest.TestCase):
    def load(self, name):
        return rt.read_fasta(FIXTURE_DIR / name)[0][1]

    def test_doubled_contig_collapses_without_anchor(self):
        seq = self.load("INT-ANIMAL-01.doubled.fasta")
        meta = json.loads(
            (FIXTURE_DIR / "INT-ANIMAL-01.doubled.bin_metadata.json")
            .read_text())
        meta = {"contigs": [
            dict(c, contig_id=CONTIG) for c in meta["contigs"]
            if c["contig_id"] == "contig_1"]}
        out, rec = self.run_tandem(seq, anchors=[], meta=meta)
        tandem = rec["tandem"]
        self.assertEqual(tandem["status"], "collapsed")
        self.assertEqual(tandem["copies"], 2)
        self.assertAlmostEqual(tandem["period_bp"], 16951, delta=5)
        self.assertGreater(tandem["support"], 0.7)
        self.assertEqual(tandem["cut_method"], "leading_kmer")
        self.assertEqual(out, seq[:tandem["cut_bp"]])

    def test_doubled_contig_rotates_to_anchor_then_cuts(self):
        seq = self.load("INT-ANIMAL-01.doubled.fasta")
        out, rec = self.run_tandem(
            seq, [hit(start=20001, end=20300)])
        self.assertEqual(rec["tandem"]["status"], "collapsed")
        self.assertAlmostEqual(len(out), 16951, delta=5)
        self.assertTrue(out.startswith(seq[20000:20031]))

    def test_single_contig_is_unchanged(self):
        seq = self.load("INT-ANIMAL-01.single.fasta")
        out, rec = self.run_tandem(seq, anchors=[])
        self.assertEqual(rec["tandem"]["status"], "single")
        self.assertLess(rec["tandem"]["support"], 0.03)
        self.assertEqual(out, seq)


class TestMain(unittest.TestCase):
    def run_main(self, tmp, anchors, gff_text):
        tmp = Path(tmp)
        seq = make_seq()
        (tmp / "target.fasta").write_text(f"{HEADER}\n{seq}\n")
        (tmp / "bin.json").write_text(json.dumps(circular_meta()))
        (tmp / "anchors.5.gff").write_text(gff_text)
        argv = [
            "rotate_target.py",
            "--target-fasta", str(tmp / "target.fasta"),
            "--bin-metadata", str(tmp / "bin.json"),
            "--anchors", anchors,
            "--gff", f"5:{tmp / 'anchors.5.gff'}",
            "--out-fasta", str(tmp / "out.fasta"),
            "--out-json", str(tmp / "rotation.json"),
        ]
        with mock.patch.object(sys, "argv", argv):
            code = rt.main()
        return code, seq, tmp

    def test_main_rotates_and_writes_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, seq, tmp = self.run_main(
                tmp, "COX1,COX2", gff_line("COX1", 101, 400) + "\n")
            self.assertEqual(code, 0)
            records = rt.read_fasta(tmp / "out.fasta")
            self.assertEqual(records[0][0], HEADER)
            self.assertEqual(records[0][1], seq[100:] + seq[:100])
            data = json.loads((tmp / "rotation.json").read_text())
            self.assertEqual(data["$schema"], "wf5/rotation/v2")
            self.assertEqual(data["status"], "rotated")
            self.assertEqual(data["anchors"], ["COX1", "COX2"])
            self.assertEqual(
                set(data["contigs"][0]),
                {"contig_id", "status", "anchor_gene", "fallback",
                 "anchor_table", "reverse_complemented", "raw_start",
                 "query_start", "offset", "length_bp", "input_length_bp",
                 "tandem", "anchors"})
            self.assertEqual(
                data["contigs"][0]["tandem"]["status"], "not_tested")

    def test_main_tandem_args_collapse(self):
        genome = make_seq(GENOME_LEN, seed=7)
        seq = genome + genome
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "target.fasta").write_text(f"{HEADER}\n{seq}\n")
            (tmp / "bin.json").write_text(json.dumps(circular_meta()))
            argv = [
                "rotate_target.py",
                "--target-fasta", str(tmp / "target.fasta"),
                "--bin-metadata", str(tmp / "bin.json"),
                "--tandem-kmer", "31",
                "--tandem-min-support", "0.5",
                "--tandem-tolerance-bp", "50",
                "--tandem-min-period-bp", "1000",
                "--out-fasta", str(tmp / "out.fasta"),
                "--out-json", str(tmp / "rotation.json"),
            ]
            with mock.patch.object(sys, "argv", argv):
                self.assertEqual(rt.main(), 0)
            data = json.loads((tmp / "rotation.json").read_text())
            contig = data["contigs"][0]
            self.assertEqual(data["tandem_params"], TANDEM)
            self.assertEqual(contig["tandem"]["status"], "collapsed")
            self.assertEqual(contig["input_length_bp"], 2 * GENOME_LEN)
            self.assertEqual(contig["length_bp"], GENOME_LEN)
            self.assertEqual(
                len(rt.read_fasta(tmp / "out.fasta")[0][1]), GENOME_LEN)

    def test_main_without_anchors(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, seq, tmp = self.run_main(tmp, "", "")
            self.assertEqual(code, 0)
            data = json.loads((tmp / "rotation.json").read_text())
            self.assertEqual(data["status"], "not_configured")
            self.assertEqual(rt.read_fasta(tmp / "out.fasta")[0][1], seq)

    def test_main_unreadable_input_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            argv = [
                "rotate_target.py",
                "--target-fasta", str(Path(tmp) / "missing.fa"),
                "--bin-metadata", str(Path(tmp) / "missing.json"),
                "--out-fasta", str(Path(tmp) / "o.fa"),
                "--out-json", str(Path(tmp) / "o.json"),
            ]
            with mock.patch.object(sys, "argv", argv):
                with contextlib.suppress(SystemExit):
                    with self.assertRaises(FileNotFoundError):
                        rt.main()


if __name__ == "__main__":
    unittest.main()
