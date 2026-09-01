"""Unit tests for bin/translate_annotation_cds.py — task 40 §2/§6.

Cases:
  1. Miniprot + rescued features both translate through the same path,
     landing in the same per-gene queries/<GENE>.faa.
  2. Minus-strand feature — reverse-complemented before translation.
  3. genetic_code_cds missing (None) — no queries/ files written
     (nothing to translate uniformly under, task 40 §2).
  4. Empty cds.gff / no CDS features parsed — no queries/ files.
  5. Feature whose spliced CDS is empty (bad coordinates) — silently
     produces no query, not dropped from the *GFF* (that never
     changes here — only the query set is affected).
  6. Two loci for the same gene (plant IR duplication) — both
     records land in that gene's single query file, each keyed by
     its own feature id.
  7. Trailing stop codon stripped from the translated protein.
"""

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN_DIR))

_spec = importlib.util.spec_from_file_location(
    "translate_annotation_cds", BIN_DIR / "translate_annotation_cds.py")
tac = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tac)

CLEAN_CODON = "GGT"  # Gly, universal under every table this repo uses


def clean_cds(n_codons: int) -> str:
    return "ATG" + CLEAN_CODON * (n_codons - 1) + "TAA"


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


class BaseCase(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _write(self, name, text):
        path = self.tmp / name
        path.write_text(text)
        return path

    def _run(self, gff_lines, fasta_records, genetic_code_cds=5):
        gff_path = self._write("merged.gff", "\n".join(gff_lines) + "\n")
        fasta_text = "".join(
            f">{seqid}\n{seq}\n" for seqid, seq in fasta_records.items()
        )
        fasta_path = self._write("target.fasta", fasta_text)
        summary_path = self._write("annotation_summary.json", json.dumps(
            {"genetic_code_cds": genetic_code_cds}))
        out_dir = self.tmp / "queries"
        tac.run(gff_path, fasta_path, summary_path, out_dir)
        return out_dir


class TranslateAnnotationCdsTests(BaseCase):

    def test_miniprot_and_rescued_share_translation_path(self):
        seq = clean_cds(10)
        out_dir = self._run(
            [
                "ctg1\tminiprot\tmRNA\t1\t{}\t.\t+\t.\t"
                "ID=MP1;Target=panel_ATP6 1 30".format(len(seq)),
                "ctg1\tminiprot\tCDS\t1\t{}\t.\t+\t0\tParent=MP1"
                .format(len(seq)),
                "ctg1\tmitos\tmRNA\t{}\t{}\t.\t+\t.\t"
                "ID=r1;Name=ATP8".format(len(seq) + 1, len(seq) * 2),
                "ctg1\tmitos\tCDS\t{}\t{}\t.\t+\t0\tParent=r1"
                .format(len(seq) + 1, len(seq) * 2),
            ],
            {"ctg1": seq + seq},
        )
        atp6 = (out_dir / "ATP6.faa").read_text()
        atp8 = (out_dir / "ATP8.faa").read_text()
        self.assertTrue(atp6.startswith(">MP1\n"))
        self.assertTrue(atp8.startswith(">r1\n"))
        # Both a miniprot winner and a rescued feature translate to
        # the same protein via extract_cds_seq + Seq.translate.
        self.assertEqual(
            atp6.splitlines()[1], atp8.splitlines()[1])

    def test_minus_strand_reverse_complemented(self):
        seq = clean_cds(5)
        rc = revcomp(seq)
        out_dir = self._run(
            [
                "ctg1\tminiprot\tmRNA\t1\t{}\t.\t-\t.\t"
                "ID=MP1;Target=panel_COX1 1 30".format(len(seq)),
                "ctg1\tminiprot\tCDS\t1\t{}\t.\t-\t0\tParent=MP1"
                .format(len(seq)),
            ],
            {"ctg1": rc},
        )
        forward_only = self._run(
            [
                "ctg1\tminiprot\tmRNA\t1\t{}\t.\t+\t.\t"
                "ID=MP1;Target=panel_COX1 1 30".format(len(seq)),
                "ctg1\tminiprot\tCDS\t1\t{}\t.\t+\t0\tParent=MP1"
                .format(len(seq)),
            ],
            {"ctg1": seq},
        )
        self.assertEqual(
            (out_dir / "COX1.faa").read_text(),
            (forward_only / "COX1.faa").read_text(),
        )

    def test_no_genetic_code_no_queries_written(self):
        seq = clean_cds(5)
        out_dir = self._run(
            [
                "ctg1\tminiprot\tmRNA\t1\t{}\t.\t+\t.\t"
                "ID=MP1;Target=panel_ATP6 1 30".format(len(seq)),
                "ctg1\tminiprot\tCDS\t1\t{}\t.\t+\t0\tParent=MP1"
                .format(len(seq)),
            ],
            {"ctg1": seq},
            genetic_code_cds=None,
        )
        self.assertEqual(list(out_dir.iterdir()), [])

    def test_no_features_no_queries_written(self):
        out_dir = self._run([], {})
        self.assertEqual(list(out_dir.iterdir()), [])

    def test_empty_splice_produces_no_query(self):
        # start > end for both mRNA and CDS is invalid but never
        # rejected upstream; extract_cds_seq slices to an empty
        # fragment and translate_feature must not choke on it.
        out_dir = self._run(
            [
                "ctg1\tminiprot\tmRNA\t5\t1\t.\t+\t.\t"
                "ID=MP1;Target=panel_ATP6 1 30",
                "ctg1\tminiprot\tCDS\t5\t1\t.\t+\t0\tParent=MP1",
            ],
            {"ctg1": "ACGT"},
        )
        self.assertEqual(list(out_dir.iterdir()), [])

    def test_two_loci_same_gene_both_in_one_file(self):
        seq = clean_cds(4)
        out_dir = self._run(
            [
                "ctg1\tminiprot\tmRNA\t1\t{}\t.\t+\t.\t"
                "ID=MP1;Target=panel_MATR 1 30".format(len(seq)),
                "ctg1\tminiprot\tCDS\t1\t{}\t.\t+\t0\tParent=MP1"
                .format(len(seq)),
                "ctg1\tminiprot\tmRNA\t{}\t{}\t.\t+\t.\t"
                "ID=MP2;Target=panel_MATR 1 30"
                .format(len(seq) + 50, len(seq) * 2 + 50),
                "ctg1\tminiprot\tCDS\t{}\t{}\t.\t+\t0\tParent=MP2"
                .format(len(seq) + 50, len(seq) * 2 + 50),
            ],
            {"ctg1": seq + ("N" * 50) + seq},
        )
        text = (out_dir / "MATR.faa").read_text()
        self.assertIn(">MP1", text)
        self.assertIn(">MP2", text)

    def test_trailing_stop_stripped(self):
        seq = clean_cds(3)
        self.assertTrue(seq.endswith("TAA"))
        out_dir = self._run(
            [
                "ctg1\tminiprot\tmRNA\t1\t{}\t.\t+\t.\t"
                "ID=MP1;Target=panel_ATP6 1 30".format(len(seq)),
                "ctg1\tminiprot\tCDS\t1\t{}\t.\t+\t0\tParent=MP1"
                .format(len(seq)),
            ],
            {"ctg1": seq},
        )
        protein = (out_dir / "ATP6.faa").read_text().splitlines()[1]
        self.assertNotIn("*", protein)

    def test_main_happy_path(self):
        seq = clean_cds(5)
        gff_path = self._write(
            "merged.gff",
            "\n".join([
                "ctg1\tminiprot\tmRNA\t1\t{}\t.\t+\t.\t"
                "ID=MP1;Target=panel_ATP6 1 30".format(len(seq)),
                "ctg1\tminiprot\tCDS\t1\t{}\t.\t+\t0\tParent=MP1"
                .format(len(seq)),
            ]) + "\n",
        )
        fasta_path = self._write("target.fasta", f">ctg1\n{seq}\n")
        summary_path = self._write(
            "annotation_summary.json",
            json.dumps({"genetic_code_cds": 5}))
        out_dir = self.tmp / "queries"

        old_argv = sys.argv
        sys.argv = [
            "translate_annotation_cds.py",
            "--merged-gff", str(gff_path),
            "--target-fasta", str(fasta_path),
            "--annotation-summary", str(summary_path),
            "--out-dir", str(out_dir),
        ]
        try:
            rc = tac.main()
        finally:
            sys.argv = old_argv
        self.assertEqual(rc, 0)
        self.assertTrue((out_dir / "ATP6.faa").exists())


if __name__ == "__main__":
    unittest.main()
