"""Unit tests for bin/validate_barcodes.py (C5) — spec §2 stage 13,
§2.2 C5, task 30 §5.

Fixtures are synthetic GFF3 + FASTA text built inline; no checked-in
binary/real miniprot output. Sequences use only codons that translate
identically under NCBI tables 1/2/5/11 (avoiding TGA/AGA/AGR) unless a
case is specifically exercising a table-dependent stop.

Cases (task 30 §5, numbered to match):
  1. Happy path — panel subset selected from a larger cds.gff.
  2. Coherence invariant — emitted coordinates match cds.gff exactly.
  3. Case-insensitive symbol match.
  4. Strand handling — minus-strand reverse complement.
  5. animal_mt clade trial — table-5-only vs valid-under-both.
  6. Neither table valid → fail, explicit reason.
  7. Internal stop codon → fail.
  8. Below identity floor → fail.
  9. Panel locus absent from cds.gff → not_found.
 10. All loci fail → empty FASTA, populated TSV, exit 0.
 11. Empty cds.gff → all loci not_found, exit 0.
 12. Duplicate gene at two loci (plastid IR case) — both emitted.
 13. Malformed GFF line → skipped with a warning, rest processed.
 14. Exit code 0 in every case — parametrised via main().

Plus: parser edge cases (attrs without '=', CDS with unknown Parent,
mRNA with no CDS children) and cluster-merge branches, for 100%
branch coverage per CONSTITUTION.md rule 14.
"""

import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[2] / "bin" / "validate_barcodes.py"
)
_spec = importlib.util.spec_from_file_location(
    "validate_barcodes", MODULE_PATH
)
vb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vb)


# Codons restricted to residues that translate identically under NCBI
# tables 1 (standard), 2 (vertebrate mt), 5 (invertebrate mt) and 11
# (plastid/bacterial) — no TGA/AGA/AGG/ATA, so a "clean" ORF is clean
# under every table this module trials.
CLEAN_CODON = "GGT"  # Gly, universal


def clean_cds(n_codons: int) -> str:
    return "ATG" + CLEAN_CODON * (n_codons - 1)


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def mrna_line(mrna_id, seqid, start, end, strand, protein_id,
              identity, extra=""):
    attrs = (
        f"ID={mrna_id};Rank=1;Identity={identity:.4f};"
        f"Target={protein_id} 1 100{extra}"
    )
    return (
        f"{seqid}\tminiprot\tmRNA\t{start}\t{end}\t100\t{strand}\t.\t"
        f"{attrs}"
    )


def cds_line(mrna_id, seqid, start, end, strand, protein_id,
             identity, extra=""):
    attrs = (
        f"Parent={mrna_id};Rank=1;Identity={identity:.4f};"
        f"Target={protein_id} 1 100{extra}"
    )
    return (
        f"{seqid}\tminiprot\tCDS\t{start}\t{end}\t100\t{strand}\t0\t"
        f"{attrs}"
    )


def feature_lines(
    mrna_id, seqid, start, end, strand, protein_id, identity, extra=""
):
    return [
        mrna_line(mrna_id, seqid, start, end, strand, protein_id,
                  identity, extra),
        cds_line(mrna_id, seqid, start, end, strand, protein_id,
                 identity, extra),
    ]


class BaseCase(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _write(self, name, text):
        path = self.tmp / name
        path.write_text(text)
        return path

    def _run(
        self, gff_lines, fasta_records, loci, assembly_target, tables,
        min_identity=60.0,
    ):
        gff_path = self._write("cds.gff", "\n".join(gff_lines) + "\n")
        fasta_text = "".join(
            f">{seqid}\n{seq}\n" for seqid, seq in fasta_records.items()
        )
        fasta_path = self._write("target.fasta", fasta_text)
        loci_path = self._write(
            "loci.json", json.dumps({assembly_target: loci}))
        out_fasta = self.tmp / "barcodes.fasta"
        out_coords = self.tmp / "coords.gff"
        out_tsv = self.tmp / "validation.tsv"

        vb.run(
            gff_path, fasta_path, assembly_target, loci_path, tables,
            min_identity, out_fasta, out_coords, out_tsv,
        )
        return out_fasta, out_coords, out_tsv

    def _tsv_rows(self, out_tsv):
        import csv
        with open(out_tsv) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def _fasta_headers(self, out_fasta):
        return [
            line[1:] for line in out_fasta.read_text().splitlines()
            if line.startswith(">")
        ]


class TestCaseMatrix(BaseCase):

    def test_case1_happy_path_subset(self):
        seq = clean_cds(30)
        # 6 panel loci, 14 off-panel genes -> 20 features total.
        panel_genes = ["COX1", "COX2", "COX3", "CYTB", "ND1", "ATP6"]
        off_panel = [f"ORF{i}" for i in range(14)]
        gff_lines = []
        seqid = "contig_1"
        genome_full = ""
        pos = 1
        for i, gene in enumerate(panel_genes + off_panel):
            start = pos
            end = start + len(seq) - 1
            genome_full += seq
            pos = end + 1
            gff_lines += feature_lines(
                f"MP{i:06d}", seqid, start, end, "+",
                f"ACC{i}_{gene}", 0.95)

        out_fasta, out_coords, out_tsv = self._run(
            gff_lines, {seqid: genome_full}, panel_genes, "animal_mt",
            [2, 5])

        headers = self._fasta_headers(out_fasta)
        self.assertEqual(len(headers), 6)
        for gene in panel_genes:
            self.assertTrue(any(h.startswith(gene + "_") for h in headers))

    def test_case2_coherence_invariant(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        start, end = 11, 10 + len(seq)
        genome = "N" * 10 + seq
        gff_lines = feature_lines(
            "MP000001", seqid, start, end, "+", "ACC1_COX1", 0.9)

        out_fasta, out_coords, out_tsv = self._run(
            gff_lines, {seqid: genome}, ["COX1"], "animal_mt", [2, 5])

        coords_text = out_coords.read_text()
        self.assertIn(
            f"{seqid}\tminiprot\tCDS\t{start}\t{end}", coords_text)
        rows = self._tsv_rows(out_tsv)
        self.assertEqual(rows[0]["seqid"], seqid)
        self.assertEqual(int(rows[0]["start"]), start)
        self.assertEqual(int(rows[0]["end"]), end)

    def test_case3_case_insensitive_symbol_match(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        gff_lines = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_cox1", 0.9)

        out_fasta, _coords, tsv = self._run(
            gff_lines, {seqid: seq}, ["COX1"], "animal_mt", [2, 5])

        self.assertEqual(self._fasta_headers(out_fasta)[0][:4], "COX1")
        self.assertEqual(self._tsv_rows(tsv)[0]["status"], "pass")

    def test_case4_strand_handling_reverse_complement(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        genome = revcomp(seq)
        gff_lines = feature_lines(
            "MP000001", seqid, 1, len(seq), "-", "ACC1_COX1", 0.9)

        out_fasta, _coords, _tsv = self._run(
            gff_lines, {seqid: genome}, ["COX1"], "animal_mt", [2, 5])

        emitted = out_fasta.read_text().splitlines()[1]
        self.assertEqual(emitted, seq)

    def test_case5a_clade_trial_table5_only(self):
        # AGA is a stop under table 2 (vertebrate) but Ser under table
        # 5 (invertebrate): valid only under 5.
        seq = "ATG" + CLEAN_CODON * 5 + "AGA" + CLEAN_CODON * 5
        seqid = "contig_1"
        gff_lines = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_COX1", 0.9)

        _fasta, _coords, tsv = self._run(
            gff_lines, {seqid: seq}, ["COX1"], "animal_mt", [2, 5])

        row = self._tsv_rows(tsv)[0]
        self.assertEqual(row["status"], "pass")
        self.assertEqual(row["genetic_code"], "5")

    def test_case5b_clade_trial_valid_under_both_picks_first(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        gff_lines = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_COX1", 0.9)

        _fasta, _coords, tsv = self._run(
            gff_lines, {seqid: seq}, ["COX1"], "animal_mt", [2, 5])

        row = self._tsv_rows(tsv)[0]
        self.assertEqual(row["genetic_code"], "2")

    def test_case6_neither_table_valid(self):
        # TAA is a universal stop under every table trialled here.
        seq = "ATG" + CLEAN_CODON * 5 + "TAA" + CLEAN_CODON * 5
        seqid = "contig_1"
        gff_lines = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_COX1", 0.9)

        _fasta, _coords, tsv = self._run(
            gff_lines, {seqid: seq}, ["COX1"], "animal_mt", [2, 5])

        row = self._tsv_rows(tsv)[0]
        self.assertEqual(row["status"], "fail")
        self.assertEqual(row["reason"], vb.REASON_INTERNAL_STOP)

    def test_case7_internal_stop_single_table(self):
        seq = "ATG" + CLEAN_CODON * 5 + "TAG" + CLEAN_CODON * 5
        seqid = "contig_1"
        gff_lines = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_rbcL", 0.9)

        _fasta, _coords, tsv = self._run(
            gff_lines, {seqid: seq}, ["rbcL"], "plant_pt", [11])

        row = self._tsv_rows(tsv)[0]
        self.assertEqual(row["status"], "fail")
        self.assertEqual(row["reason"], vb.REASON_INTERNAL_STOP)

    def test_case8_below_identity_floor(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        gff_lines = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_COX1", 0.30)

        _fasta, _coords, tsv = self._run(
            gff_lines, {seqid: seq}, ["COX1"], "animal_mt", [2, 5],
            min_identity=60.0)

        row = self._tsv_rows(tsv)[0]
        self.assertEqual(row["status"], "fail")
        self.assertEqual(row["reason"], vb.REASON_IDENTITY_BELOW_FLOOR)

    def test_case9_panel_locus_absent_not_found(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        gff_lines = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_COX1", 0.9)

        out_fasta, _coords, tsv = self._run(
            gff_lines, {seqid: seq}, ["COX1", "CYTB"], "animal_mt",
            [2, 5])

        rows = {r["gene"]: r for r in self._tsv_rows(tsv)}
        self.assertEqual(rows["CYTB"]["status"], "not_found")
        self.assertEqual(rows["CYTB"]["reason"], vb.REASON_NOT_FOUND)
        self.assertNotIn("CYTB", "".join(self._fasta_headers(out_fasta)))

    def test_case10_all_loci_fail(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        gff_lines = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_COX1", 0.1)

        out_fasta, _coords, tsv = self._run(
            gff_lines, {seqid: seq}, ["COX1"], "animal_mt", [2, 5])

        self.assertEqual(out_fasta.read_text(), "")
        rows = self._tsv_rows(tsv)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["status"], "fail")

    def test_case11_empty_cds_gff(self):
        out_fasta, out_coords, tsv = self._run(
            [], {}, ["COX1", "CYTB"], "animal_mt", [2, 5])

        self.assertEqual(out_fasta.read_text(), "")
        self.assertEqual(out_coords.read_text(), "")
        rows = self._tsv_rows(tsv)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(r["status"] == "not_found" for r in rows))

    def test_case12_duplicate_gene_two_loci_not_deduplicated(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        # Two non-overlapping loci, same gene symbol (plastid IR case).
        gff_lines = (
            feature_lines(
                "MP000001", seqid, 1, len(seq), "+", "ACC1_rbcL", 0.9)
            + feature_lines(
                "MP000002", seqid, 5000, 5000 + len(seq) - 1, "+",
                "ACC2_rbcL", 0.9)
        )
        genome = list("N" * (5000 - 1 + len(seq)))
        genome[0:len(seq)] = list(seq)
        genome[5000 - 1:5000 - 1 + len(seq)] = list(seq)
        genome_seq = "".join(genome)

        out_fasta, _coords, tsv = self._run(
            gff_lines, {seqid: genome_seq}, ["rbcL"], "plant_pt", [11])

        headers = self._fasta_headers(out_fasta)
        self.assertEqual(len(headers), 2)
        self.assertEqual(len(set(headers)), 2)
        rows = [r for r in self._tsv_rows(tsv) if r["status"] == "pass"]
        self.assertEqual(len(rows), 2)

    def test_case13_malformed_gff_line_skipped_with_warning(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        good = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_COX1", 0.9)
        gff_lines = [good[0], "not\ta\tvalid\tgff\tline"] + [good[1]]

        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            out_fasta, _coords, tsv = self._run(
                gff_lines, {seqid: seq}, ["COX1"], "animal_mt", [2, 5])

        self.assertIn("WARNING", stderr.getvalue())
        self.assertEqual(self._tsv_rows(tsv)[0]["status"], "pass")

    def test_case14_exit_code_zero_every_case(self):
        seq = clean_cds(20)
        seqid = "contig_1"
        clean_feature = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_COX1", 0.9)
        low_id_feature = feature_lines(
            "MP000001", seqid, 1, len(seq), "+", "ACC1_COX1", 0.1)
        scenarios = [
            # (gff lines, fasta, loci)
            ([], {}, ["COX1"]),
            (clean_feature, {seqid: seq}, ["COX1"]),
            (low_id_feature, {seqid: seq}, ["COX1"]),
        ]
        for gff_lines, fasta, loci in scenarios:
            with self.subTest(gff_lines=bool(gff_lines)):
                gff_path = self._write(
                    f"cds_{id(gff_lines)}.gff",
                    "\n".join(gff_lines) + "\n")
                fasta_text = "".join(
                    f">{sid}\n{s}\n" for sid, s in fasta.items())
                fasta_path = self._write(
                    f"target_{id(gff_lines)}.fasta", fasta_text)
                loci_path = self._write(
                    f"loci_{id(gff_lines)}.json",
                    json.dumps({"animal_mt": loci}))
                out_fasta = self.tmp / f"bc_{id(gff_lines)}.fasta"
                out_coords = self.tmp / f"co_{id(gff_lines)}.gff"
                out_tsv = self.tmp / f"va_{id(gff_lines)}.tsv"

                argv = [
                    "validate_barcodes.py",
                    "--cds-gff", str(gff_path),
                    "--target-fasta", str(fasta_path),
                    "--assembly-target", "animal_mt",
                    "--locus-panel", str(loci_path),
                    "--genetic-codes", "2,5",
                    "--min-identity", "60",
                    "--out-fasta", str(out_fasta),
                    "--out-coords", str(out_coords),
                    "--out-tsv", str(out_tsv),
                ]
                old_argv = sys.argv
                sys.argv = argv
                try:
                    rc = vb.main()
                finally:
                    sys.argv = old_argv
                self.assertEqual(rc, 0)


class TestParserEdgeCases(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _write(self, name, text):
        path = self.tmp / name
        path.write_text(text)
        return path

    def test_parse_attrs_skips_empty_and_keyless_tokens(self):
        attrs = vb.parse_attrs("ID=MP1;;Target=ACC_GENE 1 10;flagonly")
        self.assertEqual(attrs, {"ID": "MP1", "Target": "ACC_GENE 1 10"})

    def test_cds_with_unknown_parent_ignored(self):
        gff = self._write(
            "cds.gff",
            "contig_1\tminiprot\tCDS\t1\t30\t100\t+\t0\t"
            "Parent=NOSUCH;Target=ACC1_COX1 1 10\n",
        )
        records = vb.parse_gff(gff)
        self.assertEqual(records, [])

    def test_mrna_with_no_cds_children_filtered_out(self):
        gff = self._write(
            "cds.gff",
            mrna_line("MP000001", "contig_1", 1, 30, "+",
                      "ACC1_COX1", 0.9) + "\n",
        )
        records = vb.parse_gff(gff)
        self.assertEqual(records, [])

    def test_malformed_start_coordinate_skipped(self):
        gff = self._write(
            "cds.gff",
            "contig_1\tminiprot\tmRNA\tNOTANINT\t30\t100\t+\t.\t"
            "ID=MP1;Identity=0.9;Target=ACC1_COX1 1 10\n",
        )
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            records = vb.parse_gff(gff)
        self.assertEqual(records, [])
        self.assertIn("WARNING", stderr.getvalue())

    def test_cluster_by_locus_merges_overlapping_reps(self):
        # Two representative hits for the same gene at overlapping
        # coordinates (redundant panel representatives) collapse to
        # one winner; a third, non-overlapping hit stays separate.
        records = [
            {"gene": "COX1", "seqid": "c1", "cds": [(100, 200)],
             "identity": 0.5},
            {"gene": "COX1", "seqid": "c1", "cds": [(105, 205)],
             "identity": 0.9},
            {"gene": "COX1", "seqid": "c1", "cds": [(9000, 9100)],
             "identity": 0.4},
        ]
        winners = vb.cluster_by_locus(records)
        self.assertEqual(len(winners["COX1"]), 2)
        identities = sorted(r["identity"] for r in winners["COX1"])
        self.assertEqual(identities, [0.4, 0.9])

    def test_non_mrna_non_cds_feature_type_ignored(self):
        seq = clean_cds(20)
        lines = feature_lines(
            "MP000001", "contig_1", 1, len(seq), "+", "ACC1_COX1", 0.9)
        lines.append(
            "contig_1\tminiprot\tstop_codon\t21\t23\t0\t+\t0\t"
            "Parent=MP000001")
        gff = self._write("cds.gff", "\n".join(lines) + "\n")

        records = vb.parse_gff(gff)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["cds"], [(1, len(seq))])

    def test_extract_cds_seq_multi_exon_concatenates(self):
        from Bio.Seq import Seq

        sequences = {"c1": Seq("AAACCCGGGTTT")}
        feature = {
            "seqid": "c1", "strand": "+", "cds": [(1, 3), (7, 9)],
        }
        result = vb.extract_cds_seq(sequences, feature)
        self.assertEqual(str(result), "AAAGGG")

    def test_validate_orf_empty_sequence_invalid_length(self):
        from Bio.Seq import Seq

        passed, reason, table = vb.validate_orf(
            Seq(""), None, [1], 60.0, 0.9)
        self.assertFalse(passed)
        self.assertEqual(reason, vb.REASON_INVALID_LENGTH)
        self.assertIsNone(table)


def paf_comment(protein_id, cigar, seqid="contig_1", strand="+"):
    """Minimal ##PAF comment line carrying a cg:Z: tag — only the
    fields cigar_tag()/parse_gff() actually read need to be realistic.
    """
    return (
        f"##PAF\t{protein_id}\t100\t0\t100\t{strand}\t{seqid}\t1000\t"
        f"0\t100\t100\t100\t0\tAS:i:100\tcg:Z:{cigar}"
    )


class TestCigarAwareTranslation(BaseCase):
    """A single-base indel is the dominant real-world failure mode on
    ONT-derived assemblies (task 30 outcomes): miniprot flags it with
    a Frameshift= attribute and a cg:Z op that can't be assigned to a
    codon (F/G), rather than corrupting the reading frame outright.
    These tests pin the CIGAR-segmented translation that recovers a
    locus despite that, without ever altering the emitted sequence.
    """

    def test_cigar_tag_extracts_cg_field(self):
        line = paf_comment("ACC1_COX1", "5M2G5M")
        self.assertEqual(vb.cigar_tag(line), "5M2G5M")

    def test_cigar_tag_returns_none_when_absent(self):
        line = "##PAF\tACC1_COX1\t100\t0\t100\t+\tc1\t1000\t0\t100"
        self.assertIsNone(vb.cigar_tag(line))

    def test_parse_gff_captures_paf_cigar_for_next_mrna(self):
        seq = clean_cds(10)
        lines = [paf_comment("ACC1_COX1", "10M")] + feature_lines(
            "MP000001", "contig_1", 1, len(seq), "+", "ACC1_COX1", 0.9)
        gff = self._write("cds.gff", "\n".join(lines) + "\n")

        records = vb.parse_gff(gff)

        self.assertEqual(records[0]["cigar"], "10M")

    def test_parse_gff_cigar_not_carried_to_unrelated_mrna(self):
        seq = clean_cds(10)
        lines = (
            [paf_comment("ACC1_COX1", "10M")]
            + feature_lines(
                "MP000001", "contig_1", 1, len(seq), "+",
                "ACC1_COX1", 0.9)
            + feature_lines(
                "MP000002", "contig_1", 5000, 5000 + len(seq) - 1,
                "+", "ACC2_COX1", 0.9)
        )
        gff = self._write("cds.gff", "\n".join(lines) + "\n")

        records = {r["protein_id"]: r for r in vb.parse_gff(gff)}

        self.assertEqual(records["ACC1_COX1"]["cigar"], "10M")
        self.assertIsNone(records["ACC2_COX1"]["cigar"])

    def test_codon_blocks_frame_safe_ops_only(self):
        # 2M (6nt) + 1D (3nt, in-frame deletion, still real codons) +
        # 1I (0nt, protein-only) + 2M (6nt) — all frame-safe, one block.
        seq = "AAACCCGGGTTTGGGCCC"
        blocks = vb.codon_blocks(seq, "2M1D1I2M")
        self.assertEqual(blocks, [seq])

    def test_codon_blocks_frameshift_break_drops_bases(self):
        # A single-base deletion error: 5 clean codons, a 2-nt
        # frameshift junction that can't form a codon, 5 more clean
        # codons. The 2 junction bases are excluded from both blocks.
        block1 = clean_cds(5)
        block2 = clean_cds(5)
        garbage = "TA"
        genome_seq = block1 + garbage + block2

        blocks = vb.codon_blocks(genome_seq, "5M2G5M")

        self.assertEqual(blocks, [block1, block2])

    def test_codon_blocks_intron_ops_break_like_frameshift(self):
        block1 = clean_cds(4)
        block2 = clean_cds(4)
        intron = "GTAAGT"  # 6nt, arbitrary — dropped either way
        genome_seq = block1 + intron + block2

        blocks = vb.codon_blocks(genome_seq, "4M6N4M")

        self.assertEqual(blocks, [block1, block2])

    def test_codon_blocks_trailing_codon_remainder_kept(self):
        # A trailing stop codon appended past the CIGAR's own span
        # (task 30 outcomes: cg:Z covers the protein-coding portion,
        # not the terminal stop) is folded into the last block.
        body = clean_cds(4)
        stop = "TAA"
        seq = body + stop

        blocks = vb.codon_blocks(seq, "4M")

        self.assertEqual(blocks, [body + stop])

    def test_codon_blocks_non_codon_remainder_dropped(self):
        body = clean_cds(4)
        seq = body + "TA"  # 2 leftover nt — not a whole codon

        blocks = vb.codon_blocks(seq, "4M")

        self.assertEqual(blocks, [body])

    def test_block_has_internal_stop_last_block_strips_trailing(self):
        codon_seq = clean_cds(3) + "TAA"
        self.assertFalse(
            vb.block_has_internal_stop(codon_seq, 1, is_last=True))

    def test_block_has_internal_stop_middle_block_keeps_trailing(self):
        # A stop at the end of a non-final block is not a trailing
        # terminator — it's a real internal stop.
        codon_seq = clean_cds(3) + "TAA"
        self.assertTrue(
            vb.block_has_internal_stop(codon_seq, 1, is_last=False))

    def test_validate_orf_recovers_single_base_indel_via_cigar(self):
        # Naive whole-sequence translation hits a spurious stop at the
        # frameshift junction; CIGAR-segmented translation does not.
        block1 = clean_cds(5)
        block2 = clean_cds(5)
        genome_seq = block1 + "TA" + block2  # "TA"+"A.." -> TAA naively

        naive = vb.validate_orf(genome_seq, None, [1], 60.0, 0.95)
        cigar_aware = vb.validate_orf(
            genome_seq, "5M2G5M", [1], 60.0, 0.95)

        self.assertFalse(naive[0])
        self.assertEqual(naive[1], vb.REASON_INTERNAL_STOP)
        self.assertTrue(cigar_aware[0])

    def test_validate_orf_cigar_with_no_frame_safe_blocks_fails(self):
        # Every op is a block-break op — nothing left to translate.
        passed, reason, table = vb.validate_orf(
            "TAGTC", "5G", [1], 60.0, 0.95)
        self.assertFalse(passed)
        self.assertEqual(reason, vb.REASON_INVALID_LENGTH)

    def test_end_to_end_single_base_indel_recovered_unmodified_seq(self):
        # Full run() pipeline: emitted barcode sequence is the raw
        # coordinate-verbatim genome slice, unmodified by the CIGAR —
        # only the pass/fail verdict is CIGAR-aware.
        block1 = clean_cds(5)
        block2 = clean_cds(5)
        genome_seq = block1 + "TA" + block2
        seqid = "contig_1"
        lines = [paf_comment("ACC1_COX1", "5M2G5M")] + feature_lines(
            "MP000001", seqid, 1, len(genome_seq), "+", "ACC1_COX1",
            0.9)

        out_fasta, _coords, tsv = self._run(
            lines, {seqid: genome_seq}, ["COX1"], "animal_mt", [1])

        self.assertEqual(self._tsv_rows(tsv)[0]["status"], "pass")
        emitted = out_fasta.read_text().splitlines()[1]
        self.assertEqual(emitted, genome_seq)


if __name__ == "__main__":
    unittest.main()
