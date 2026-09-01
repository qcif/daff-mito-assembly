"""Unit tests for bin/select_genetic_code.py — spec §2 stage 12,
task 38 §2/§3.

Fixtures are small hand-written miniprot-shaped GFF3 fragments, not
captured tool output (CONSTITUTION rule 19).

Cases:
  1. Single-candidate passthrough — one table, no real decision.
  2. Two-candidate selection — table A wins on total score.
  3. Two-candidate selection — table B wins (opposite direction).
  4. Exact tie (equal score, equal gene count) — earliest-listed wins.
  5. Empty candidate GFF (miniprot found nothing under one table).
  6. Malformed input: wrong column count, non-numeric score, missing
     Target attribute — each skipped with a warning, selection still
     correct.
  7. main() end-to-end via parse_candidate_arg + CLI.

Plus direct branch coverage of gene_from_target / parse_attrs /
parse_best_scores per CONSTITUTION rule 14.
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
    Path(__file__).resolve().parents[2] / "bin" / "select_genetic_code.py"
)
_spec = importlib.util.spec_from_file_location(
    "select_genetic_code", MODULE_PATH
)
sgc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sgc)


def mrna_line(seqid, start, end, score, gene, strand="+", n=1):
    target = f"NC_000000_{gene} 1 100"
    return (
        f"{seqid}\tminiprot\tmRNA\t{start}\t{end}\t{score}\t{strand}\t.\t"
        f"ID=MP{n:06d};Identity=1.0;Target={target}"
    )


def write_gff(tmp: Path, name: str, lines: list) -> Path:
    path = tmp / name
    path.write_text("\n".join(lines) + "\n" if lines else "")
    return path


class SelectGeneticCodeTests(unittest.TestCase):

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmpdir.name)

    def tearDown(self):
        self._tmpdir.cleanup()

    def run_select(self, candidates):
        """candidates: [(table, [lines])] -> (winner, records)."""
        pairs = []
        for i, (table, lines) in enumerate(candidates):
            path = write_gff(self.tmp, f"cds.{table}.{i}.gff", lines)
            pairs.append((table, path))
        return sgc.select(pairs)

    # -- 1. single-candidate passthrough --------------------------------

    def test_single_candidate_passthrough(self):
        winner, records = self.run_select([
            (11, [mrna_line("ctg1", 1, 300, 500, "RBCL")]),
        ])
        self.assertEqual(len(records), 1)
        self.assertEqual(winner["table"], 11)
        self.assertEqual(winner["total_score"], 500.0)
        self.assertEqual(winner["n_genes"], 1)

    # -- 2/3. two-candidate selection, either direction ------------------

    def test_table_a_wins_on_score(self):
        winner, records = self.run_select([
            (2, [mrna_line("ctg1", 1, 300, 300, "COX1")]),
            (5, [mrna_line("ctg1", 1, 300, 700, "COX1")]),
        ])
        self.assertEqual(winner["table"], 5)
        self.assertEqual(len(records), 2)

    def test_table_b_wins_on_score(self):
        winner, _ = self.run_select([
            (2, [mrna_line("ctg1", 1, 300, 900, "COX1")]),
            (5, [mrna_line("ctg1", 1, 300, 100, "COX1")]),
        ])
        self.assertEqual(winner["table"], 2)

    def test_tie_break_on_gene_count(self):
        winner, _ = self.run_select([
            (2, [mrna_line("ctg1", 1, 300, 500, "COX1")]),
            (5, [
                mrna_line("ctg1", 1, 300, 300, "COX1", n=1),
                mrna_line("ctg1", 400, 600, 200, "COX2", n=2),
            ]),
        ])
        # Equal total_score (500 == 300+200); table 5 wins on n_genes.
        self.assertEqual(winner["table"], 5)

    # -- 4. exact tie -----------------------------------------------------

    def test_exact_tie_keeps_earliest_listed(self):
        winner, records = self.run_select([
            (2, [mrna_line("ctg1", 1, 300, 500, "COX1")]),
            (5, [mrna_line("ctg1", 1, 300, 500, "COX1")]),
        ])
        self.assertEqual(winner["table"], 2)
        self.assertEqual(records[0]["total_score"], records[1]["total_score"])
        self.assertEqual(records[0]["n_genes"], records[1]["n_genes"])

    # -- 5. empty candidate GFF --------------------------------------------

    def test_empty_candidate_loses(self):
        winner, records = self.run_select([
            (2, []),
            (5, [mrna_line("ctg1", 1, 300, 500, "COX1")]),
        ])
        self.assertEqual(winner["table"], 5)
        empty_record = next(r for r in records if r["table"] == 2)
        self.assertEqual(empty_record["total_score"], 0)
        self.assertEqual(empty_record["n_genes"], 0)

    def test_both_candidates_empty(self):
        winner, records = self.run_select([(2, []), (5, [])])
        self.assertEqual(winner["table"], 2)
        self.assertEqual(winner["total_score"], 0)
        self.assertEqual(winner["n_genes"], 0)

    # -- 6. malformed input -------------------------------------------------

    def test_wrong_column_count_skipped(self):
        stderr = io.StringIO()
        path = write_gff(self.tmp, "cds.bad.gff", [
            "ctg1\tminiprot\tmRNA\t1\t300",  # too few columns
            mrna_line("ctg1", 1, 300, 400, "COX1"),
        ])
        with contextlib.redirect_stderr(stderr):
            scores = sgc.parse_best_scores(path)
        self.assertEqual(scores, {"COX1": 400.0})

    def test_non_numeric_score_skipped_with_warning(self):
        stderr = io.StringIO()
        path = write_gff(self.tmp, "cds.badscore.gff", [
            "ctg1\tminiprot\tmRNA\t1\t300\tNOTANUM\t+\t.\t"
            "ID=MP1;Target=NC_1_COX1 1 100",
            mrna_line("ctg1", 1, 300, 400, "COX2"),
        ])
        with contextlib.redirect_stderr(stderr):
            scores = sgc.parse_best_scores(path)
        self.assertIn("non-numeric score", stderr.getvalue())
        self.assertEqual(scores, {"COX2": 400.0})

    def test_missing_target_skipped_with_warning(self):
        stderr = io.StringIO()
        path = write_gff(self.tmp, "cds.notarget.gff", [
            "ctg1\tminiprot\tmRNA\t1\t300\t400\t+\t.\tID=MP1",
            mrna_line("ctg1", 1, 300, 400, "COX2"),
        ])
        with contextlib.redirect_stderr(stderr):
            scores = sgc.parse_best_scores(path)
        self.assertIn("no Target attribute", stderr.getvalue())
        self.assertEqual(scores, {"COX2": 400.0})

    def test_non_mrna_rows_ignored(self):
        path = write_gff(self.tmp, "cds.withcds.gff", [
            mrna_line("ctg1", 1, 300, 400, "COX1"),
            "ctg1\tminiprot\tCDS\t1\t300\t400\t+\t0\tParent=MP1",
        ])
        self.assertEqual(sgc.parse_best_scores(path), {"COX1": 400.0})

    def test_blank_and_comment_lines_ignored(self):
        path = write_gff(self.tmp, "cds.comment.gff", [
            "##gff-version 3",
            "",
            mrna_line("ctg1", 1, 300, 400, "COX1"),
        ])
        self.assertEqual(sgc.parse_best_scores(path), {"COX1": 400.0})

    def test_second_lower_scoring_hit_for_same_gene_not_kept(self):
        path = write_gff(self.tmp, "cds.dup.gff", [
            mrna_line("ctg1", 1, 300, 900, "COX1", n=1),
            mrna_line("ctg1", 400, 700, 300, "COX1", n=2),
        ])
        self.assertEqual(sgc.parse_best_scores(path), {"COX1": 900.0})

    def test_second_higher_scoring_hit_for_same_gene_replaces_first(self):
        path = write_gff(self.tmp, "cds.dup2.gff", [
            mrna_line("ctg1", 1, 300, 300, "COX1", n=1),
            mrna_line("ctg1", 400, 700, 900, "COX1", n=2),
        ])
        self.assertEqual(sgc.parse_best_scores(path), {"COX1": 900.0})

    # -- helpers ------------------------------------------------------------

    def test_gene_from_target(self):
        self.assertEqual(
            sgc.gene_from_target("NC_012345_COX1 1 100"), "COX1")

    def test_parse_attrs(self):
        self.assertEqual(
            sgc.parse_attrs("ID=MP1;Target=NC_1_COX1 1 100"),
            {"ID": "MP1", "Target": "NC_1_COX1 1 100"},
        )

    def test_parse_attrs_ignores_malformed_kv(self):
        self.assertEqual(
            sgc.parse_attrs("ID=MP1;novalue;Target=NC_1_COX1 1 100"),
            {"ID": "MP1", "Target": "NC_1_COX1 1 100"},
        )

    # -- 7. main() / CLI ------------------------------------------------------

    def test_parse_candidate_arg(self):
        table, path = sgc.parse_candidate_arg("5:some/dir/cds.5.gff")
        self.assertEqual(table, 5)
        self.assertEqual(path, Path("some/dir/cds.5.gff"))

    def test_main_end_to_end(self):
        gff_a = write_gff(self.tmp, "cds.2.gff", [
            mrna_line("ctg1", 1, 300, 300, "COX1"),
        ])
        gff_b = write_gff(self.tmp, "cds.5.gff", [
            mrna_line("ctg1", 1, 300, 700, "COX1"),
        ])
        out_gff = self.tmp / "winner.gff"
        out_json = self.tmp / "genetic_code.json"

        argv = [
            "select_genetic_code.py",
            "--candidate", f"2:{gff_a}",
            "--candidate", f"5:{gff_b}",
            "--out-gff", str(out_gff),
            "--out-json", str(out_json),
        ]
        old_argv = sys.argv
        sys.argv = argv
        try:
            rc = sgc.main()
        finally:
            sys.argv = old_argv

        self.assertEqual(rc, 0)
        self.assertEqual(out_gff.read_text(), gff_b.read_text())
        result = json.loads(out_json.read_text())
        self.assertEqual(result["selected_table"], 5)
        self.assertEqual(len(result["candidates"]), 2)
        self.assertIn("criterion", result)


if __name__ == "__main__":
    unittest.main()
