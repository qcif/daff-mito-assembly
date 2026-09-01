"""Unit tests for bin/annotate_summary.py (C8) — spec §2 stage 13a,
§2.2 C8, task 31 §5.

Fixtures are small on-disk GFF/JSON trees, not mocks — C8 shells out
to nothing. A synthetic gene-sets JSON (not the real
assets/organelle_gene_sets.json) keeps these tests independent of
reference-data drift (CONSTITUTION rule 19).

Cases (task 31 §5, numbered to match):
  1. Happy path (animal_mt) — cds.gff + a MITOS2 tree merge.
  2. CDS provenance — every CDS feature traces to a cds.gff row.
  3. ok_cds_only — no annotator directory supplied.
  4. no_assembly — empty cds.gff.
  5. no_features — both sources present but nothing parseable.
  6. Cross-check agreement.
  7. Cross-check disagreement (single-source + coordinate conflict).
  8. IR duplication retained.
  9. Name normalisation.
 10. Multi-record input (MITOS2 per-record subdirectories).
 11. Genetic-code mismatch.
 12. Malformed GFF line skipped with a warning.
 13. Unknown assembly_target.
 14. Exit code 0 in every case.
 15. annotator_failed — configured annotator, no result.gff at all
     (task 41 §5.1).
 16. annotator_failed — result.gff present but parses to zero
     tRNA/rRNA rows (task 41 §5.2).
 17. Not confused with ok_cds_only when no annotator is configured
     at all (task 41 §5.3).
 18. annotator_exit_code recorded when given, null when omitted
     (task 41 §5.4).
 19. ok is unchanged on the task 31 happy-path fixture (task 41 §5.5).
 20. A legitimate rRNA-only (no tRNA) annotation is still ok, not
     annotator_failed — the trigger is zero features, not a missing
     type (task 41 §5.6).
 21. Rescue of annotator-only CDS calls and its four guards
     (task 39 §2/§6) — see TestRescueAnnotatorOnly.

Plus parser/helper edge-case branches for 100% branch coverage per
CONSTITUTION rule 14.
"""

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[2] / "bin" / "annotate_summary.py"
)
_spec = importlib.util.spec_from_file_location(
    "annotate_summary", MODULE_PATH
)
asum = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(asum)


GENE_SETS = {
    "$schema": "wf5/gene-sets/v1",
    "version": "test",
    "sets": {
        "animal_mt": {
            "name": "metazoan_37",
            "protein_coding": [
                "COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8",
            ],
            "protein_coding_aliases": {"CYTB": ["COB", "CYB"]},
            "rrna": ["rrnL", "rrnS"],
            "trna_count": 22,
        },
        "empty_target": {
            "name": "empty",
            "protein_coding": [],
            "protein_coding_aliases": {},
            "rrna": [],
        },
    },
}


def cds_gff_text(rows):
    """rows: list of (gene, seqid, start, end, strand, identity)."""
    lines = ["##gff-version 3"]
    for i, (gene, seqid, start, end, strand, identity) in enumerate(
            rows, 1):
        target = f"NC_000001.1_{gene} 1 100"
        lines.append(
            f"{seqid}\tminiprot\tmRNA\t{start}\t{end}\t{identity}\t"
            f"{strand}\t.\tID=MP{i:06d};Identity={identity};"
            f"Target={target}"
        )
        lines.append(
            f"{seqid}\tminiprot\tCDS\t{start}\t{end}\t{identity}\t"
            f"{strand}\t0\tParent=MP{i:06d};Identity={identity};"
            f"Target={target}"
        )
    return "\n".join(lines) + "\n"


def mitos_gff_text(seqid, rows):
    """rows: list of (ftype, name, start, end, strand) for tRNA/rRNA/
    gene rows, or (ftype, name, start, end, strand, extra) for an
    ``exon`` row — ``extra`` is a dict of ``phase`` (default ``"0"``),
    ``parent`` (default ``transcript_<name>``, matching MITOS2's own
    naming convention) and ``exon_name`` (default ``name``)."""
    lines = ["##gff-version 3", "#!gff-spec-version 1.21"]
    for row in rows:
        ftype, name, start, end, strand = row[:5]
        extra = row[5] if len(row) > 5 else {}
        if ftype == "gene":
            lines.append(
                f"{seqid}\tmitos\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
                f"ID=gene_{name};Name={name};gene_id={name}"
            )
        elif ftype == "exon":
            parent = extra.get("parent", f"transcript_{name}")
            phase = extra.get("phase", "0")
            exon_name = extra.get("exon_name", name)
            lines.append(
                f"{seqid}\tmitos\texon\t{start}\t{end}\t.\t{strand}\t"
                f"{phase}\tID=exon_{exon_name};Parent={parent};"
                f"Name={exon_name}"
            )
        else:
            lines.append(
                f"{seqid}\tmitfi\t{ftype}\t{start}\t{end}\t.\t{strand}\t.\t"
                f"ID=transcript_{name};Name={name};Parent=gene_{name};"
                f"gene_id={name}"
            )
    return "\n".join(lines) + "\n"


class RunHelper(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.gene_sets_path = self.dir / "gene_sets.json"
        self.gene_sets_path.write_text(json.dumps(GENE_SETS))

    def write_cds_gff(self, rows, name="cds.gff"):
        path = self.dir / name
        path.write_text(cds_gff_text(rows))
        return path

    def write_mitos_dir(self, seqid, rows, subdir=None):
        mitos_dir = self.dir / "mitos_out"
        target = mitos_dir / subdir if subdir else mitos_dir
        target.mkdir(parents=True, exist_ok=True)
        (target / "result.gff").write_text(mitos_gff_text(seqid, rows))
        return mitos_dir

    def run_c8(self, cds_gff, mitos_dir=None, assembly_target="animal_mt",
               genetic_code_annotate=None, genetic_code_cds=None,
               reference_data=None, annotator_exit=None,
               sample_id="SAMPLE-01"):
        out_gff = self.dir / "out.gff"
        out_summary = self.dir / "out.summary.json"
        asum.run(
            cds_gff, sample_id, assembly_target, self.gene_sets_path,
            mitos_dir, genetic_code_annotate, genetic_code_cds,
            reference_data, annotator_exit, out_gff, out_summary,
        )
        summary = json.loads(out_summary.read_text())
        gff_text = out_gff.read_text()
        return summary, gff_text


class TestHappyPath(RunHelper):
    """Case 1 — cds.gff + a MITOS2 tree merge, status ok."""

    def test_happy_path(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
            ("ATP6", "contig_1", 500, 700, "+", 0.85),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("tRNA", "trnF(gaa)", 10, 80, "+"),
            ("rRNA", "rrnS", 800, 1200, "+"),
            ("gene", "cox1", 100, 400, "+"),
        ])
        summary, gff_text = self.run_c8(
            cds_gff, mitos_dir=mitos_dir, genetic_code_annotate=5,
            genetic_code_cds=5, reference_data="refseq89m")

        self.assertEqual(summary["status"], "ok")
        self.assertEqual(summary["reason"], None)
        self.assertEqual(
            summary["cds_source"], {"miniprot": 2, "mitos": 0})
        self.assertEqual(summary["non_cds_source"], "mitos2")
        self.assertEqual(summary["reference_data"], "refseq89m")
        self.assertEqual(
            summary["feature_counts"],
            {"gene": 4, "CDS": 2, "tRNA": 1, "rRNA": 1})
        self.assertEqual(summary["genetic_code_agreement"], True)
        self.assertIn("contig_1", summary["contigs_annotated"])
        self.assertIn("trnF(gaa)", gff_text)
        self.assertIn("rrnS", gff_text)


class TestCdsProvenance(RunHelper):
    """Case 2 — every CDS feature traces to a cds.gff row exactly; no
    MITOS2 CDS ('gene' type) rows ever appear as features."""

    def test_no_mitos_cds_leaks_into_gff(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "nad2", 900, 1000, "+"),
        ])
        summary, gff_text = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertNotIn("nad2", gff_text)
        self.assertIn("100\t400", gff_text)
        self.assertEqual(summary["feature_counts"]["CDS"], 1)


class TestOkCdsOnly(RunHelper):
    """Case 3 — no annotator directory supplied."""

    def test_ok_cds_only(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=None)
        self.assertEqual(summary["status"], "ok_cds_only")
        self.assertIsNotNone(summary["reason"])
        self.assertIsNone(summary["non_cds_source"])
        self.assertIsNone(summary["cds_crosscheck"])
        self.assertIsNone(summary["cds_rescued"])


class TestNoAssembly(RunHelper):
    """Case 4 — empty cds.gff, exit 0."""

    def test_empty_cds_gff(self):
        cds_gff = self.dir / "empty.gff"
        cds_gff.write_text("")
        summary, gff_text = self.run_c8(cds_gff)
        self.assertEqual(summary["status"], "no_assembly")
        self.assertEqual(summary["feature_counts"]["CDS"], 0)
        self.assertEqual(gff_text.strip(), "##gff-version 3")

    def test_missing_cds_gff(self):
        cds_gff = self.dir / "does_not_exist.gff"
        summary, _ = self.run_c8(cds_gff)
        self.assertEqual(summary["status"], "no_assembly")


class TestNoFeatures(RunHelper):
    """Case 5 — non-empty cds.gff but nothing parseable; no annotator
    features either -> no_features."""

    def test_no_features_parsed(self):
        cds_gff = self.dir / "junk.gff"
        cds_gff.write_text("not\ta\tvalid\tgff\tline\n")
        summary, _ = self.run_c8(cds_gff)
        self.assertEqual(summary["status"], "no_features")


class TestCrossCheck(RunHelper):
    """Cases 6/7 — agreement, single-source, coordinate conflict."""

    def test_agreement(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "cox1", 150, 380, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertIn("cox1", summary["cds_crosscheck"]["agreed"])

    def test_miniprot_only_and_annotator_only(self):
        # nad2 is deliberately absent from this trimmed test gene set,
        # so it also exercises the §2.3 off-panel guard.
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "nad2", 900, 1000, "+"),
            ("exon", "nad2", 900, 1000, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        cc = summary["cds_crosscheck"]
        self.assertIn("COX1", cc["miniprot_only"])
        self.assertEqual(
            cc["annotator_only"],
            [{"gene": "nad2", "reason": "off_panel"}])
        self.assertEqual(cc["agreed"], [])
        self.assertEqual(cc["coordinate_conflicts"], [])
        self.assertEqual(summary["cds_rescued"], [])

    def test_coordinate_conflict(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "cox1", 5000, 5300, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        cc = summary["cds_crosscheck"]
        self.assertIn("cox1", cc["coordinate_conflicts"])
        self.assertEqual(cc["agreed"], [])
        self.assertEqual(cc["miniprot_only"], [])
        self.assertEqual(cc["annotator_only"], [])


class TestRescueAnnotatorOnly(RunHelper):
    """Task 39 — an ``annotator_only`` protein-coding call with no
    overlapping miniprot CDS is emitted as a real feature, attributed
    to its annotator, and counted toward completeness."""

    def test_clean_single_exon_rescue(self):
        """The ATP8 case (task 39 §0/§1): a short, fast-evolving gene
        miniprot never saw at all, MITOS2 called correctly."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp8", 1719, 1919, "+"),
            ("exon", "atp8", 1719, 1919, "+", {"phase": "0"}),
        ])
        summary, gff_text = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        cc = summary["cds_crosscheck"]
        self.assertEqual(cc["annotator_only"], [])
        self.assertEqual(summary["cds_rescued"], ["ATP8"])
        self.assertNotIn("ATP8", summary["protein_coding_genes_missing"])
        self.assertEqual(summary["feature_counts"]["CDS"], 2)
        self.assertEqual(
            summary["cds_source"], {"miniprot": 1, "mitos": 1})
        self.assertIn("contig_1\tmitos\tmRNA\t1719\t1919", gff_text)
        self.assertIn("contig_1\tmitos\tCDS\t1719\t1919", gff_text)

    def test_multi_exon_rescue_preserves_phase(self):
        """§2.2 — MITOS2 puts the gene span on the ``gene`` row and
        the reading frame on ``exon`` rows; rescue must assemble the
        CDS from the exons, not the gene row."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp6", 1197, 1615, "+"),
            ("exon", "atp6", 1197, 1481, "+",
             {"phase": "0", "exon_name": "atp6-b"}),
            ("exon", "atp6", 1496, 1615, "+",
             {"phase": "2", "exon_name": "atp6-a"}),
        ])
        summary, gff_text = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(summary["cds_rescued"], ["ATP6"])
        self.assertIn("contig_1\tmitos\tCDS\t1197\t1481\t.\t+\t0", gff_text)
        self.assertIn("contig_1\tmitos\tCDS\t1496\t1615\t.\t+\t2", gff_text)

    def test_fragment_suffix_folds_and_does_not_duplicate(self):
        """§2.1 — overlapping fragment-suffixed calls (``atp8_0``/
        ``atp8_1``) fold to one canonical name and rescue as a single
        feature, not two."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp8_1", 1719, 1919, "+"),
            ("exon", "atp8_1", 1719, 1919, "+",
             {"phase": "0", "parent": "transcript_atp8_1"}),
            ("gene", "atp8_0", 1897, 2031, "+"),
            ("exon", "atp8_0", 1897, 2031, "+",
             {"phase": "0", "parent": "transcript_atp8_0"}),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(summary["cds_rescued"], ["ATP8"])
        self.assertEqual(summary["feature_counts"]["CDS"], 2)

    def test_off_panel_call_held_back(self):
        """§2.3 — a call outside the canonical set (e.g. a LAGLIDADG
        homing endonuclease) is never rescued."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "lagli", 5000, 5300, "+"),
            ("exon", "lagli", 5000, 5300, "+"),
        ])
        summary, gff_text = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        cc = summary["cds_crosscheck"]
        self.assertEqual(
            cc["annotator_only"],
            [{"gene": "lagli", "reason": "off_panel"}])
        self.assertEqual(summary["cds_rescued"], [])
        self.assertNotIn("lagli", gff_text)

    def test_overlapping_call_held_back(self):
        """§2.4 — a coordinate guard, not just a name test: a call
        that folds to a name with no miniprot counterpart but sits on
        top of an existing miniprot CDS on the same strand is held
        back, not rescued as a duplicate-by-position feature."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp8", 150, 380, "+"),
            ("exon", "atp8", 150, 380, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        cc = summary["cds_crosscheck"]
        self.assertEqual(
            cc["annotator_only"],
            [{"gene": "atp8", "reason": "overlap"}])
        self.assertEqual(summary["cds_rescued"], [])

    def test_off_panel_beats_overlap_as_the_reported_reason(self):
        """A call that is both off-panel and positionally colliding
        reports off_panel — the disqualifier that would hold it back
        whatever its coordinates. INT-ANIMAL-01's lagli sits wholly
        inside miniprot's COX1; reporting 'overlap' would point an
        auditor at the wrong problem."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "lagli", 150, 300, "+"),
            ("exon", "lagli", 150, 300, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(
            summary["cds_crosscheck"]["annotator_only"],
            [{"gene": "lagli", "reason": "off_panel"}])

    def test_small_junction_overlap_still_rescues(self):
        """§2.4 — a few bp of overlap with a neighbouring miniprot CDS
        is normal in a compact organelle genome (§1 cites the canonical
        7 bp ATP8/ATP6 junction as evidence the gene is *real*), so it
        must not be mistaken for duplication. 7 bp of a 207 bp call is
        3.4% of the shorter feature, far below the threshold."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp8", 394, 600, "+"),
            ("exon", "atp8", 394, 600, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(summary["cds_crosscheck"]["annotator_only"], [])
        self.assertEqual(summary["cds_rescued"], ["ATP8"])

    def test_overlap_threshold_boundary(self):
        """Pins RESCUE_MAX_OVERLAP_FRACTION itself: 50/100 of the
        shorter feature is rejected (>=), 49/100 rescues."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        at_threshold = self.write_mitos_dir("contig_1", [
            ("gene", "atp8", 351, 450, "+"),
            ("exon", "atp8", 351, 450, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=at_threshold)
        self.assertEqual(summary["cds_rescued"], [])
        self.assertEqual(
            summary["cds_crosscheck"]["annotator_only"],
            [{"gene": "atp8", "reason": "overlap"}])

        below = self.write_mitos_dir("contig_1", [
            ("gene", "atp8", 352, 451, "+"),
            ("exon", "atp8", 352, 451, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=below)
        self.assertEqual(summary["cds_rescued"], ["ATP8"])

    def test_colliding_fragment_does_not_suppress_clean_sibling(self):
        """§2.4 guard order, from the real INT-ANIMAL-01 geometry:
        MITOS2 offers atp8_1 (overlapping miniprot's over-extended
        ATP6 by 64%) and atp8_0 (clean). Clustering before the guard
        picks atp8_1 on length and discards atp8_0, losing ATP8
        entirely; guarding first keeps the rescuable one."""
        cds_gff = self.write_cds_gff([
            ("ATP6", "contig_1", 1200, 1847, "-", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp8_1", 1719, 1919, "-"),
            ("exon", "atp8_1", 1719, 1919, "-",
             {"parent": "transcript_atp8_1"}),
            ("gene", "atp8_0", 1897, 2031, "-"),
            ("exon", "atp8_0", 1897, 2031, "-",
             {"parent": "transcript_atp8_0"}),
        ])
        summary, gff_text = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(summary["cds_rescued"], ["ATP8"])
        self.assertEqual(
            summary["cds_crosscheck"]["annotator_only"],
            [{"gene": "atp8_1", "reason": "overlap"}])
        self.assertIn("contig_1\tmitos\tCDS\t1897\t2031", gff_text)
        self.assertNotIn("1719\t1919", gff_text)

    def test_different_overlapping_genes_both_survive(self):
        """§2.1 — clustering is per canonical gene. Two *different*
        genes abutting on the same strand are two genes; a
        position-only cluster would drop the shorter with no feature
        and no annotator_only record at all (principle 18)."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp8", 5000, 5200, "+"),
            ("exon", "atp8", 5000, 5200, "+"),
            ("gene", "atp6", 5150, 5600, "+"),
            ("exon", "atp6", 5150, 5600, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(summary["cds_rescued"], ["ATP6", "ATP8"])
        self.assertEqual(summary["cds_crosscheck"]["annotator_only"], [])

    def test_overlap_guard_is_strand_specific(self):
        """Same coordinates as the held-back case above, but on the
        opposite strand from the miniprot CDS — dense organelle
        genomes legitimately carry genes on complementary strands at
        overlapping coordinates, so this must rescue, not hold back."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp8", 150, 380, "-"),
            ("exon", "atp8", 150, 380, "-"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(summary["cds_crosscheck"]["annotator_only"], [])
        self.assertEqual(summary["cds_rescued"], ["ATP8"])

    def test_no_exon_data_held_back(self):
        """Boundary guard: a gene row with no matching exon row (not
        MITOS2's documented shape, but this process doesn't control
        MITOS2's output) is held back rather than emitted with no
        CDS lines."""
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp8", 1719, 1919, "+"),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        cc = summary["cds_crosscheck"]
        self.assertEqual(
            cc["annotator_only"],
            [{"gene": "atp8", "reason": "no_exon_data"}])
        self.assertEqual(summary["cds_rescued"], [])

    def test_completeness_before_and_after_rescue(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        no_mitos_summary, _ = self.run_c8(cds_gff, mitos_dir=None)
        self.assertIn(
            "ATP8", no_mitos_summary["protein_coding_genes_missing"])

        mitos_dir = self.write_mitos_dir("contig_1", [
            ("gene", "atp8", 1719, 1919, "+"),
            ("exon", "atp8", 1719, 1919, "+"),
        ])
        rescued_summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertNotIn(
            "ATP8", rescued_summary["protein_coding_genes_missing"])
        self.assertGreater(
            rescued_summary["protein_coding_completeness"],
            no_mitos_summary["protein_coding_completeness"])


class TestIRDuplication(RunHelper):
    """Case 8 — same gene at two non-overlapping ranges survives as
    two winners, both retained in the GFF, counted as two loci."""

    def test_duplicate_gene_two_loci(self):
        cds_gff = self.write_cds_gff([
            ("psbA", "path1", 100, 400, "+", 0.9),
            ("psbA", "path1", 9000, 9300, "+", 0.9),
        ])
        summary, gff_text = self.run_c8(
            cds_gff, assembly_target="animal_mt")
        self.assertEqual(summary["feature_counts"]["CDS"], 2)
        self.assertEqual(gff_text.count("psbA"), 4)  # 2x mRNA + 2x CDS

    def test_overlapping_hits_collapse_to_one_winner(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.8),
            ("COX1", "contig_1", 110, 410, "+", 0.95),
        ])
        summary, gff_text = self.run_c8(cds_gff)
        self.assertEqual(summary["feature_counts"]["CDS"], 1)
        self.assertIn("110\t410", gff_text)
        self.assertNotIn("100\t400", gff_text)


class TestNameNormalisation(RunHelper):
    """Case 9 — NAD4L / trnF(gaa) / rrnS style names match the
    canonical set for completeness/cross-check purposes; raw GFF
    attribute text is left untouched."""

    def test_alias_and_anticodon_folding(self):
        cds_gff = self.write_cds_gff([
            ("CYTB", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("tRNA", "trnF(gaa)", 10, 80, "+"),
            ("gene", "cob", 100, 400, "+"),
        ])
        summary, gff_text = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertIn("cob", summary["cds_crosscheck"]["agreed"])
        self.assertIn("trnF(gaa)", gff_text)
        self.assertNotIn("CYTB", summary["protein_coding_genes_missing"])

    def test_mt_prefix_stripped(self):
        self.assertEqual(
            asum.normalise_gene_symbol("MT-CO1"), "CO1")

    def test_fragment_suffix_stripped(self):
        self.assertEqual(asum.normalise_gene_symbol("atp8_0"), "ATP8")
        self.assertEqual(asum.normalise_gene_symbol("cox1_1"), "COX1")

    def test_fragment_suffix_does_not_eat_real_trailing_digits(self):
        self.assertEqual(asum.normalise_gene_symbol("nad4l"), "NAD4L")


class TestMultiRecord(RunHelper):
    """Case 10 — plant_mt's emit: all case; MITOS2 subdirectories
    (one per input FASTA record) unioned, each already carrying the
    correct real seqid in column 1."""

    def test_multi_record_subdirs(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_a", 100, 400, "+", 0.9),
            ("ATP6", "contig_b", 50, 200, "+", 0.9),
        ])
        mitos_dir = self.dir / "mitos_out"
        (mitos_dir / "0").mkdir(parents=True)
        (mitos_dir / "0" / "result.gff").write_text(
            mitos_gff_text("contig_a", [("tRNA", "trnF(gaa)", 10, 80, "+")])
        )
        (mitos_dir / "1").mkdir(parents=True)
        (mitos_dir / "1" / "result.gff").write_text(
            mitos_gff_text("contig_b", [("rRNA", "rrnS", 700, 900, "+")])
        )
        summary, gff_text = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(summary["feature_counts"]["tRNA"], 1)
        self.assertEqual(summary["feature_counts"]["rRNA"], 1)
        self.assertIn("contig_a\tmitfi\ttRNA", gff_text)
        self.assertIn("contig_b\tmitfi\trRNA", gff_text)
        self.assertEqual(
            set(summary["contigs_annotated"]), {"contig_a", "contig_b"})


class TestGeneticCodeAgreement(RunHelper):
    """Case 11 — MINIPROT_CDS's actual per-sample selection compared
    against MITOS2's fixed configured table (task 38 §4). Unlike the
    pre-task-38 config-constant proxy, this comparison can genuinely
    disagree: a vertebrate animal_mt submission where MINIPROT_CDS
    selects table 2 while MITOS2 stays configured at 5 is a real QC
    signal, not a bug in the comparison."""

    def test_mismatch(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        summary, _ = self.run_c8(
            cds_gff, genetic_code_annotate=5, genetic_code_cds=2)
        self.assertEqual(summary["genetic_code_agreement"], False)
        self.assertNotEqual(summary["status"], "no_assembly")

    def test_agreement(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        summary, _ = self.run_c8(
            cds_gff, genetic_code_annotate=5, genetic_code_cds=5)
        self.assertEqual(summary["genetic_code_agreement"], True)

    def test_one_side_missing_is_none(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        summary, _ = self.run_c8(
            cds_gff, genetic_code_annotate=5, genetic_code_cds=None)
        self.assertIsNone(summary["genetic_code_agreement"])


class TestMalformedGff(RunHelper):
    """Case 12 — malformed lines skipped with a warning; the rest of
    the file is still processed."""

    def test_malformed_line_skipped(self):
        cds_gff = self.dir / "cds.gff"
        good = cds_gff_text([("COX1", "contig_1", 100, 400, "+", 0.9)])
        cds_gff.write_text(good + "too\tfew\tcols\n")
        summary, _ = self.run_c8(cds_gff)
        self.assertEqual(summary["feature_counts"]["CDS"], 1)

    def test_non_integer_coordinates_skipped(self):
        cds_gff = self.dir / "cds.gff"
        cds_gff.write_text(
            "contig_1\tminiprot\tmRNA\tNOTANUM\t400\t0.9\t+\t.\t"
            "ID=MP000001;Identity=0.9;Target=NC_1_COX1 1 100\n"
        )
        summary, _ = self.run_c8(cds_gff)
        self.assertEqual(summary["status"], "no_features")

    def test_mrna_missing_target_skipped(self):
        cds_gff = self.dir / "cds.gff"
        cds_gff.write_text(
            "contig_1\tminiprot\tmRNA\t100\t400\t0.9\t+\t.\tID=MP1\n"
        )
        summary, _ = self.run_c8(cds_gff)
        self.assertEqual(summary["feature_counts"]["CDS"], 0)

    def test_cds_with_unknown_parent_skipped(self):
        cds_gff = self.dir / "cds.gff"
        cds_gff.write_text(
            "contig_1\tminiprot\tCDS\t100\t400\t0.9\t+\t0\t"
            "Parent=MPnotreal;Identity=0.9\n"
        )
        summary, _ = self.run_c8(cds_gff)
        self.assertEqual(summary["status"], "no_features")

    def test_cds_gff_ignores_other_feature_types(self):
        cds_gff = self.dir / "cds.gff"
        good = cds_gff_text([("COX1", "contig_1", 100, 400, "+", 0.9)])
        cds_gff.write_text(
            good + "contig_1\tminiprot\tstop_codon\t401\t403\t.\t+\t0\t"
            "Parent=MP000001\n"
        )
        summary, _ = self.run_c8(cds_gff)
        self.assertEqual(summary["feature_counts"]["CDS"], 1)

    def test_mitos_malformed_line_skipped(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.dir / "mitos_out"
        mitos_dir.mkdir()
        (mitos_dir / "result.gff").write_text(
            mitos_gff_text("contig_1", [("tRNA", "trnF(gaa)", 10, 80, "+")])
            + "broken\tline\n"
        )
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(summary["feature_counts"]["tRNA"], 1)

    def test_mitos_exon_without_parent_skipped(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.dir / "mitos_out"
        mitos_dir.mkdir()
        (mitos_dir / "result.gff").write_text(
            mitos_gff_text("contig_1", [("gene", "atp8", 1719, 1919, "+")])
            + "contig_1\tmitos\texon\t1719\t1919\t.\t+\t0\tID=exon_atp8\n"
        )
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        cc = summary["cds_crosscheck"]
        self.assertEqual(
            cc["annotator_only"],
            [{"gene": "atp8", "reason": "no_exon_data"}])

    def test_mitos_gff_ignores_other_feature_types(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.dir / "mitos_out"
        mitos_dir.mkdir()
        (mitos_dir / "result.gff").write_text(
            mitos_gff_text("contig_1", [("tRNA", "trnF(gaa)", 10, 80, "+")])
            + "contig_1\tmitos\tregion\t1\t16799\t.\t+\t.\tID=contig_1\n"
        )
        summary, _ = self.run_c8(cds_gff, mitos_dir=mitos_dir)
        self.assertEqual(summary["feature_counts"]["tRNA"], 1)


class TestAnnotatorFailed(RunHelper):
    """Cases 15-20 — task 41: a configured non-CDS annotator that
    yields zero features is reported, never silently folded into
    ``ok``."""

    def test_no_result_gff_at_all(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.dir / "mitos_out"
        mitos_dir.mkdir()
        summary, gff_text = self.run_c8(
            cds_gff, mitos_dir=mitos_dir, annotator_exit=2)
        self.assertEqual(summary["status"], "annotator_failed")
        self.assertIn("mitos2", summary["reason"])
        self.assertIn("2", summary["reason"])
        self.assertIn("mitos_out", summary["reason"])
        # CDS features are still emitted in full.
        self.assertEqual(summary["feature_counts"]["CDS"], 1)
        self.assertIn("100\t400", gff_text)

    def test_empty_result_gff(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [])
        summary, _ = self.run_c8(
            cds_gff, mitos_dir=mitos_dir, annotator_exit=0)
        self.assertEqual(summary["status"], "annotator_failed")
        self.assertEqual(summary["feature_counts"]["CDS"], 1)

    def test_not_confused_with_ok_cds_only(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=None)
        self.assertEqual(summary["status"], "ok_cds_only")
        self.assertIsNone(summary["non_cds_source"])

    def test_annotator_exit_code_recorded(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("tRNA", "trnF(gaa)", 10, 80, "+"),
        ])
        summary, _ = self.run_c8(
            cds_gff, mitos_dir=mitos_dir, annotator_exit=2)
        self.assertEqual(summary["annotator_exit_code"], 2)

    def test_annotator_exit_code_null_when_omitted(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        summary, _ = self.run_c8(cds_gff, mitos_dir=None)
        self.assertIsNone(summary["annotator_exit_code"])

    def test_happy_path_still_ok(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("tRNA", "trnF(gaa)", 10, 80, "+"),
            ("rRNA", "rrnS", 800, 1200, "+"),
        ])
        summary, _ = self.run_c8(
            cds_gff, mitos_dir=mitos_dir, annotator_exit=0)
        self.assertEqual(summary["status"], "ok")

    def test_rrna_only_is_ok_not_failed(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        mitos_dir = self.write_mitos_dir("contig_1", [
            ("rRNA", "rrnS", 800, 1200, "+"),
        ])
        summary, _ = self.run_c8(
            cds_gff, mitos_dir=mitos_dir, annotator_exit=0)
        self.assertEqual(summary["status"], "ok")


class TestUnknownAssemblyTarget(RunHelper):
    """Case 13 — unknown assembly_target -> completeness fields None,
    status unaffected."""

    def test_unknown_target(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        summary, _ = self.run_c8(cds_gff, assembly_target="not_a_target")
        self.assertIsNone(summary["protein_coding_completeness"])
        self.assertIsNone(summary["protein_coding_genes_missing"])
        self.assertIsNone(summary["canonical_gene_set"])
        self.assertEqual(summary["status"], "ok_cds_only")

    def test_empty_canonical_set(self):
        cds_gff = self.write_cds_gff([
            ("COX1", "contig_1", 100, 400, "+", 0.9),
        ])
        summary, _ = self.run_c8(cds_gff, assembly_target="empty_target")
        self.assertIsNone(summary["protein_coding_completeness"])


class TestExitCodeZero(unittest.TestCase):
    """Case 14 — exit code 0 in every case, parametrised via main()."""

    def _run_main(self, argv):
        old_argv = sys.argv
        sys.argv = ["annotate_summary.py"] + argv
        try:
            return asum.main()
        finally:
            sys.argv = old_argv

    def test_main_happy_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            gene_sets_path = d / "gene_sets.json"
            gene_sets_path.write_text(json.dumps(GENE_SETS))
            cds_gff = d / "cds.gff"
            cds_gff.write_text(
                cds_gff_text([("COX1", "contig_1", 100, 400, "+", 0.9)]))
            genetic_code_json = d / "genetic_code.json"
            genetic_code_json.write_text(json.dumps({"selected_table": 5}))
            out_gff = d / "out.gff"
            out_summary = d / "out.summary.json"
            rc = self._run_main([
                "--cds-gff", str(cds_gff),
                "--sample-id", "S1",
                "--assembly-target", "animal_mt",
                "--gene-sets", str(gene_sets_path),
                "--genetic-code-annotate", "5",
                "--genetic-code-json", str(genetic_code_json),
                "--reference-data", "refseq89m",
                "--out-gff", str(out_gff),
                "--out-summary", str(out_summary),
            ])
            self.assertEqual(rc, 0)
            self.assertTrue(out_summary.exists())
            summary = json.loads(out_summary.read_text())
            self.assertEqual(summary["genetic_code_cds"], 5)

    def test_main_annotator_exit_arg(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            gene_sets_path = d / "gene_sets.json"
            gene_sets_path.write_text(json.dumps(GENE_SETS))
            cds_gff = d / "cds.gff"
            cds_gff.write_text(
                cds_gff_text([("COX1", "contig_1", 100, 400, "+", 0.9)]))
            mitos_dir = d / "mitos_out"
            mitos_dir.mkdir()
            out_gff = d / "out.gff"
            out_summary = d / "out.summary.json"
            rc = self._run_main([
                "--cds-gff", str(cds_gff),
                "--sample-id", "S1",
                "--assembly-target", "animal_mt",
                "--gene-sets", str(gene_sets_path),
                "--mitos-dir", str(mitos_dir),
                "--annotator-exit", "2",
                "--out-gff", str(out_gff),
                "--out-summary", str(out_summary),
            ])
            self.assertEqual(rc, 0)
            summary = json.loads(out_summary.read_text())
            self.assertEqual(summary["status"], "annotator_failed")
            self.assertEqual(summary["annotator_exit_code"], 2)

    def test_main_no_optional_args(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            gene_sets_path = d / "gene_sets.json"
            gene_sets_path.write_text(json.dumps(GENE_SETS))
            cds_gff = d / "cds.gff"
            cds_gff.write_text("")
            out_gff = d / "out.gff"
            out_summary = d / "out.summary.json"
            rc = self._run_main([
                "--cds-gff", str(cds_gff),
                "--sample-id", "S1",
                "--assembly-target", "animal_mt",
                "--gene-sets", str(gene_sets_path),
                "--out-gff", str(out_gff),
                "--out-summary", str(out_summary),
            ])
            self.assertEqual(rc, 0)
            summary = json.loads(out_summary.read_text())
            self.assertEqual(summary["status"], "no_assembly")


class TestHelperUnits(unittest.TestCase):
    """Direct branch coverage for small helpers not otherwise exercised
    through run()."""

    def test_parse_attrs_skips_bad_pairs(self):
        self.assertEqual(
            asum.parse_attrs("a=1;;noequals;b=2"), {"a": "1", "b": "2"})

    def test_parse_gff_line_wrong_column_count(self):
        self.assertIsNone(
            asum.parse_gff_line("a\tb\tc", 1, Path("x.gff")))

    def test_overlaps_false_different_seqid(self):
        a = {"seqid": "c1", "start": 1, "end": 10}
        b = {"seqid": "c2", "start": 1, "end": 10}
        self.assertFalse(asum.overlaps(a, b))

    def test_overlaps_false_no_range_overlap(self):
        a = {"seqid": "c1", "start": 1, "end": 10}
        b = {"seqid": "c1", "start": 20, "end": 30}
        self.assertFalse(asum.overlaps(a, b))

    def test_build_alias_index_no_aliases_key(self):
        gene_set = {"protein_coding": ["COX1"]}
        index = asum.build_alias_index(gene_set)
        self.assertEqual(index, {"COX1": "COX1"})


if __name__ == "__main__":
    unittest.main()
