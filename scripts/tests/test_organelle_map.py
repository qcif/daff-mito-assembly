"""Unit tests for bin/organelle_map.py — task 44 (stage 14 ORGANELLE_MAP).

Cases:
  1. parse_non_cds_features: tRNA/rRNA parsed; Name/gene_id/type gene
     fallback chain; malformed/non-9-column/non-numeric lines skipped.
  2. load_cds_scores: missing path, zero-byte, malformed JSON, normal.
  3. cds_features: gene from cds_scores; falls back to the GFF-derived
     gene when the id has no score entry.
  4. cluster_loci: singleton non-CDS features never cluster; genomically
     overlapping same-strand CDS calls collapse to the best-scored
     primary with the rest recorded as alt_genes; different strand or
     non-overlapping calls stay separate.
  5. genome_length: plastid substitution path1_len; contigs_selected
     fallback; max-feature-end fallback; None when nothing is known.
  6. plastid_path2_segments: no canonicalisation block, substitution not
     applied, incomplete segment lengths, and the valid case.
  7. remap_to_path2: feature entirely in LSC/IR unchanged; feature in
     SSC mirrored + strand-flipped; boundary-spanning feature dropped.
  8. build_panels / run: two-panel output when path2 applies; empty
     output file when genome length cannot be determined.
  9. render_svg: wedge geometry — large-arc branch, minimum-wedge-width
     branch, +/- strand ring placement, tooltip with/without alt_genes
     and score fields.
 10. main(): end-to-end success, and the catch-all exception path still
     exits 0 and leaves an empty file.
 11. (task 54 §5.4) load_locus_panel / cluster_loci barcode-panel
     labelling: case-insensitive match, per-target panel selection,
     primary-gene-only (never alt_genes), tRNA/rRNA never labelled,
     unreadable/absent panel gives no labels, plastid path2 carries
     labels through the remap and drops boundary-spanning features.
 12. (task 54 §5.3) legend/label containment — every legend and label
     text position lies within the rendered viewBox.
"""

import importlib.util
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN_DIR))

_spec = importlib.util.spec_from_file_location(
    "organelle_map", BIN_DIR / "organelle_map.py")
om = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(om)


def _write(tmp: Path, name: str, text: str) -> Path:
    path = tmp / name
    path.write_text(text)
    return path


class ParseNonCdsFeaturesTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_trna_and_rrna_parsed(self):
        gff = _write(self.tmp, "a.gff", (
            "##gff-version 3\n"
            "contig_1\tmitfi\ttRNA\t10\t20\t.\t+\t.\t"
            "ID=t1;Name=trnG(tcc);Parent=gene_trnG\n"
            "contig_1\tmitfi\trRNA\t30\t40\t.\t-\t.\t"
            "ID=r1;Name=rrnS;Parent=gene_rrnS\n"
        ))
        feats = om.parse_non_cds_features(gff)
        self.assertEqual(len(feats), 2)
        self.assertEqual(feats[0]['kind'], 'tRNA')
        self.assertEqual(feats[0]['gene'], 'trnG(tcc)')
        self.assertEqual(feats[1]['strand'], '-')

    def test_gene_id_fallback_when_no_name(self):
        gff = _write(self.tmp, "a.gff", (
            "contig_1\tmitfi\ttRNA\t10\t20\t.\t+\t.\t"
            "ID=t1;gene_id=trnQ\n"
        ))
        feats = om.parse_non_cds_features(gff)
        self.assertEqual(feats[0]['gene'], 'trnQ')

    def test_type_fallback_when_no_name_or_gene_id(self):
        gff = _write(self.tmp, "a.gff", (
            "contig_1\tmitfi\trRNA\t10\t20\t.\t+\t.\tID=r1\n"
        ))
        feats = om.parse_non_cds_features(gff)
        self.assertEqual(feats[0]['gene'], 'rRNA')

    def test_skips_non_9_column_and_comment_and_blank_lines(self):
        gff = _write(self.tmp, "a.gff", (
            "##gff-version 3\n"
            "\n"
            "contig_1\tmitfi\ttRNA\t10\t20\t.\t+\t.\n"
            "contig_1\tmitfi\ttRNA\t10\t20\t.\t+\t.\tID=t1;Name=trnA\n"
        ))
        feats = om.parse_non_cds_features(gff)
        self.assertEqual(len(feats), 1)

    def test_skips_non_trna_rrna_types(self):
        gff = _write(self.tmp, "a.gff", (
            "contig_1\tminiprot\tmRNA\t10\t20\t.\t+\t.\tID=m1\n"
        ))
        self.assertEqual(om.parse_non_cds_features(gff), [])

    def test_skips_non_numeric_coordinates(self):
        gff = _write(self.tmp, "a.gff", (
            "contig_1\tmitfi\ttRNA\tNA\t20\t.\t+\t.\tID=t1;Name=trnA\n"
        ))
        self.assertEqual(om.parse_non_cds_features(gff), [])

    def test_attribute_without_equals_sign_ignored(self):
        gff = _write(self.tmp, "a.gff", (
            "contig_1\tmitfi\ttRNA\t10\t20\t.\t+\t.\t"
            "malformed;ID=t1;Name=trnA\n"
        ))
        feats = om.parse_non_cds_features(gff)
        self.assertEqual(feats[0]['gene'], 'trnA')


class LoadCdsScoresTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_missing_path_returns_empty(self):
        self.assertEqual(om.load_cds_scores(None), {})

    def test_missing_file_returns_empty(self):
        self.assertEqual(
            om.load_cds_scores(self.tmp / "nope.json"), {})

    def test_zero_byte_file_returns_empty(self):
        path = _write(self.tmp, "empty.json", "")
        self.assertEqual(om.load_cds_scores(path), {})

    def test_malformed_json_returns_empty(self):
        path = _write(self.tmp, "bad.json", "{not json")
        self.assertEqual(om.load_cds_scores(path), {})

    def test_normal_summary(self):
        path = _write(self.tmp, "summary.json", json.dumps({
            "cds_scores": [
                {"id": "MP1", "gene": "COX1", "pident": 90.0},
                {"id": "MP2", "gene": "ND1"},
            ]
        }))
        scores = om.load_cds_scores(path)
        self.assertEqual(scores["MP1"]["gene"], "COX1")
        self.assertEqual(scores["MP2"]["gene"], "ND1")

    def test_row_without_id_skipped(self):
        path = _write(self.tmp, "summary.json", json.dumps({
            "cds_scores": [{"gene": "COX1"}]
        }))
        self.assertEqual(om.load_cds_scores(path), {})


class LoadBinMetadataTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_none_path_returns_empty(self):
        self.assertEqual(om.load_bin_metadata(None), {})

    def test_missing_file_returns_empty(self):
        self.assertEqual(
            om.load_bin_metadata(self.tmp / "nope.json"), {})

    def test_zero_byte_file_returns_empty(self):
        path = _write(self.tmp, "empty.json", "")
        self.assertEqual(om.load_bin_metadata(path), {})

    def test_malformed_json_returns_empty(self):
        path = _write(self.tmp, "bad.json", "{not json")
        self.assertEqual(om.load_bin_metadata(path), {})

    def test_normal_metadata(self):
        path = _write(
            self.tmp, "meta.json", json.dumps({'contigs_selected': []}))
        self.assertEqual(
            om.load_bin_metadata(path), {'contigs_selected': []})


class CdsFeaturesTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_gene_from_score_overrides_gff(self):
        gff = _write(self.tmp, "a.gff", (
            "contig_1\tminiprot\tmRNA\t10\t50\t.\t+\t.\t"
            "ID=MP1;Target=NC_1_cox1 1 10\n"
            "contig_1\tminiprot\tCDS\t10\t50\t.\t+\t0\tParent=MP1\n"
        ))
        scores = {"MP1": {
            "gene": "COX1", "pident": 90.0, "qcovhsp": 80.0,
            "bitscore": 100.0}}
        feats = om.cds_features(gff, scores)
        self.assertEqual(feats[0]['gene'], 'COX1')
        self.assertEqual(feats[0]['pident'], 90.0)

    def test_falls_back_to_gff_gene_when_unscored(self):
        gff = _write(self.tmp, "a.gff", (
            "contig_1\tminiprot\tmRNA\t10\t50\t.\t+\t.\t"
            "ID=MP1;Target=NC_1_cox1 1 10\n"
            "contig_1\tminiprot\tCDS\t10\t50\t.\t+\t0\tParent=MP1\n"
        ))
        feats = om.cds_features(gff, {})
        self.assertEqual(feats[0]['gene'], 'cox1')
        self.assertIsNone(feats[0]['pident'])


class ClusterLociTests(unittest.TestCase):

    def _cds(self, gene, start, end, strand='+', pident=None,
             bitscore=None):
        return {
            'seqid': 'c1', 'source': 'miniprot', 'kind': 'CDS',
            'strand': strand, 'start': start, 'end': end, 'gene': gene,
            'pident': pident, 'qcovhsp': None, 'bitscore': bitscore,
        }

    def test_singleton_non_cds_not_clustered(self):
        other = [{
            'seqid': 'c1', 'source': 'mitfi', 'kind': 'tRNA',
            'strand': '+', 'start': 1, 'end': 10, 'gene': 'trnA',
            'pident': None, 'qcovhsp': None, 'bitscore': None,
        }]
        loci = om.cluster_loci([], other)
        self.assertEqual(len(loci), 1)
        self.assertEqual(loci[0]['alt_genes'], [])

    def test_overlapping_same_strand_collapse_to_best_scored(self):
        feats = [
            self._cds('COX1', 10, 100, bitscore=50.0),
            self._cds('cox1', 10, 100, pident=90.0, bitscore=200.0),
            self._cds('coI', 15, 90),
        ]
        loci = om.cluster_loci(feats, [])
        self.assertEqual(len(loci), 1)
        self.assertEqual(loci[0]['gene'], 'cox1')
        self.assertEqual(set(loci[0]['alt_genes']), {'COX1', 'coI'})

    def test_non_overlapping_stays_separate(self):
        feats = [
            self._cds('ND1', 10, 50),
            self._cds('ND2', 1000, 1050),
        ]
        loci = om.cluster_loci(feats, [])
        self.assertEqual(len(loci), 2)

    def test_different_strand_not_clustered(self):
        feats = [
            self._cds('ND1', 10, 50, strand='+'),
            self._cds('ND1', 10, 50, strand='-'),
        ]
        loci = om.cluster_loci(feats, [])
        self.assertEqual(len(loci), 2)

    def test_tie_break_prefers_longer_span_when_unscored(self):
        feats = [
            self._cds('short', 10, 20),
            self._cds('long', 10, 100),
        ]
        loci = om.cluster_loci(feats, [])
        self.assertEqual(loci[0]['gene'], 'long')

    def test_cluster_group_of_empty_list_returns_no_clusters(self):
        self.assertEqual(om._cluster_group([]), [])

    def test_other_feats_never_flagged_as_barcode(self):
        other = [{
            'seqid': 'c1', 'source': 'mitfi', 'kind': 'tRNA',
            'strand': '+', 'start': 1, 'end': 10, 'gene': 'COX1',
            'pident': None, 'qcovhsp': None, 'bitscore': None,
        }]
        loci = om.cluster_loci([], other, frozenset({'COX1'}))
        self.assertFalse(loci[0]['is_barcode'])

    def test_primary_gene_in_panel_is_labelled(self):
        feats = [self._cds('COX1', 10, 100)]
        loci = om.cluster_loci(feats, [], frozenset({'COX1'}))
        self.assertTrue(loci[0]['is_barcode'])

    def test_non_panel_gene_not_labelled(self):
        feats = [self._cds('ND1', 10, 100)]
        loci = om.cluster_loci(feats, [], frozenset({'COX1'}))
        self.assertFalse(loci[0]['is_barcode'])

    def test_case_insensitive_match(self):
        feats = [self._cds('rbcL', 10, 100)]
        loci = om.cluster_loci(feats, [], frozenset({'RBCL'}))
        self.assertTrue(loci[0]['is_barcode'])

    def test_panel_gene_only_in_alt_genes_not_labelled(self):
        # The cluster's primary (highest score) is 'ND1'; 'COX1' is
        # only an alt_gene at this locus. The panel gene appearing as
        # an alt_gene must not label the arc (§5.4 — labelled for what
        # it was drawn as).
        feats = [
            self._cds('ND1', 10, 100, bitscore=200.0),
            self._cds('COX1', 10, 100, bitscore=50.0),
        ]
        loci = om.cluster_loci(feats, [], frozenset({'COX1'}))
        self.assertEqual(len(loci), 1)
        self.assertEqual(loci[0]['gene'], 'ND1')
        self.assertIn('COX1', loci[0]['alt_genes'])
        self.assertFalse(loci[0]['is_barcode'])

    def test_empty_panel_labels_nothing(self):
        feats = [self._cds('COX1', 10, 100)]
        loci = om.cluster_loci(feats, [], frozenset())
        self.assertFalse(loci[0]['is_barcode'])


class LoadLocusPanelTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_none_path_or_target_returns_empty(self):
        self.assertEqual(om.load_locus_panel(None, 'animal_mt'), frozenset())
        path = _write(self.tmp, "p.json", json.dumps({"animal_mt": ["COX1"]}))
        self.assertEqual(om.load_locus_panel(path, None), frozenset())

    def test_missing_file_returns_empty(self):
        self.assertEqual(
            om.load_locus_panel(self.tmp / "nope.json", "animal_mt"),
            frozenset())

    def test_zero_byte_file_returns_empty(self):
        path = _write(self.tmp, "empty.json", "")
        self.assertEqual(
            om.load_locus_panel(path, "animal_mt"), frozenset())

    def test_malformed_json_returns_empty(self):
        path = _write(self.tmp, "bad.json", "{not json")
        self.assertEqual(
            om.load_locus_panel(path, "animal_mt"), frozenset())

    def test_target_with_no_panel_entry_returns_empty(self):
        path = _write(
            self.tmp, "panel.json", json.dumps({"animal_mt": ["COX1"]}))
        self.assertEqual(
            om.load_locus_panel(path, "plant_pt"), frozenset())

    def test_normal_panel_uppercased(self):
        path = _write(
            self.tmp, "panel.json",
            json.dumps({"plant_pt": ["rbcL", "matK"]}))
        self.assertEqual(
            om.load_locus_panel(path, "plant_pt"), frozenset({
                "RBCL", "MATK"}))

    def test_different_target_genes_do_not_cross_label(self):
        # An animal_mt panel gene must not label a plant_pt map.
        path = _write(
            self.tmp, "panel.json", json.dumps({
                "animal_mt": ["COX1"], "plant_pt": ["rbcL"],
            }))
        panel = om.load_locus_panel(path, "plant_pt")
        self.assertNotIn("COX1", panel)
        self.assertIn("RBCL", panel)


class GenomeLengthTests(unittest.TestCase):

    def test_substitution_applied_uses_path1_len(self):
        metadata = {
            'plastid_canonicalisation': {
                'substitution_applied': True, 'path1_len': 12345,
            },
        }
        self.assertEqual(om.genome_length(metadata, []), 12345)

    def test_contigs_selected_fallback(self):
        metadata = {
            'contigs_selected': ['contig_3'],
            'contigs': [
                {'contig_id': 'contig_1', 'length_bp': 999},
                {'contig_id': 'contig_3', 'length_bp': 16500},
            ],
        }
        self.assertEqual(om.genome_length(metadata, []), 16500)

    def test_falls_back_to_max_feature_end(self):
        loci = [{'end': 500}, {'end': 900}]
        self.assertEqual(om.genome_length({}, loci), 900)

    def test_none_when_nothing_available(self):
        self.assertIsNone(om.genome_length({}, []))

    def test_selected_contig_with_falsy_length_falls_through(self):
        metadata = {
            'contigs_selected': ['c1'],
            'contigs': [{'contig_id': 'c1', 'length_bp': 0}],
        }
        loci = [{'end': 777}]
        self.assertEqual(om.genome_length(metadata, loci), 777)

    def test_substitution_flag_without_path1_len_falls_through(self):
        metadata = {
            'plastid_canonicalisation': {'substitution_applied': True},
            'contigs_selected': ['c1'],
            'contigs': [{'contig_id': 'c1', 'length_bp': 42}],
        }
        self.assertEqual(om.genome_length(metadata, []), 42)


class PlastidPath2SegmentsTests(unittest.TestCase):

    def test_no_canonicalisation_block(self):
        self.assertIsNone(om.plastid_path2_segments({}))

    def test_substitution_not_applied(self):
        metadata = {'plastid_canonicalisation': {
            'substitution_applied': False}}
        self.assertIsNone(om.plastid_path2_segments(metadata))

    def test_incomplete_segment_lengths(self):
        metadata = {'plastid_canonicalisation': {
            'substitution_applied': True, 'lsc_len': 90000,
            'ir_len': None, 'ssc_len': 15000}}
        self.assertIsNone(om.plastid_path2_segments(metadata))

    def test_valid_segments(self):
        metadata = {'plastid_canonicalisation': {
            'substitution_applied': True, 'lsc_len': 90000,
            'ir_len': 25000, 'ssc_len': 15000}}
        self.assertEqual(
            om.plastid_path2_segments(metadata), (90000, 25000, 15000))


class RemapToPath2Tests(unittest.TestCase):

    SEGMENTS = (100, 20, 30)  # lsc=100, ir=20, ssc=30 -> ssc at [121,150]

    def _locus(self, start, end, strand='+'):
        return {
            'seqid': 'path1', 'source': 'miniprot', 'kind': 'CDS',
            'strand': strand, 'start': start, 'end': end, 'gene': 'g',
            'pident': None, 'qcovhsp': None, 'bitscore': None,
            'alt_genes': [],
        }

    def test_feature_outside_ssc_unchanged(self):
        locus = self._locus(10, 50)
        [remapped] = om.remap_to_path2([locus], self.SEGMENTS)
        self.assertEqual((remapped['start'], remapped['end']), (10, 50))
        self.assertEqual(remapped['strand'], '+')

    def test_feature_in_final_ir_unchanged(self):
        locus = self._locus(160, 165)
        [remapped] = om.remap_to_path2([locus], self.SEGMENTS)
        self.assertEqual((remapped['start'], remapped['end']), (160, 165))

    def test_feature_inside_ssc_mirrored_and_flipped(self):
        # SSC occupies 1-based [121, 150]. A feature at [121, 130]
        # (the first 10 bases of the SSC) mirrors to the last 10 bases.
        locus = self._locus(121, 130, strand='+')
        [remapped] = om.remap_to_path2([locus], self.SEGMENTS)
        self.assertEqual((remapped['start'], remapped['end']), (141, 150))
        self.assertEqual(remapped['strand'], '-')

    def test_mirror_is_involutive(self):
        locus = self._locus(121, 130)
        [once] = om.remap_to_path2([locus], self.SEGMENTS)
        [twice] = om.remap_to_path2([once], self.SEGMENTS)
        self.assertEqual((twice['start'], twice['end']), (121, 130))
        self.assertEqual(twice['strand'], '+')

    def test_feature_spanning_ssc_boundary_dropped(self):
        locus = self._locus(115, 125)  # crosses IR/SSC boundary at 120
        self.assertEqual(om.remap_to_path2([locus], self.SEGMENTS), [])

    def test_unknown_strand_left_unflipped(self):
        locus = self._locus(121, 130, strand='.')
        [remapped] = om.remap_to_path2([locus], self.SEGMENTS)
        self.assertEqual(remapped['strand'], '.')

    def test_is_barcode_flag_carried_through_remap(self):
        locus = {**self._locus(10, 50), 'is_barcode': True}
        [remapped] = om.remap_to_path2([locus], self.SEGMENTS)
        self.assertTrue(remapped['is_barcode'])

    def test_boundary_spanning_barcode_locus_dropped_from_path2(self):
        locus = {**self._locus(115, 125), 'is_barcode': True}
        self.assertEqual(om.remap_to_path2([locus], self.SEGMENTS), [])


class BuildPanelsAndRunTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def _gff(self):
        return _write(self.tmp, "a.gff", (
            "contig_1\tminiprot\tmRNA\t10\t100\t.\t+\t.\t"
            "ID=MP1;Target=NC_1_cox1 1 30\n"
            "contig_1\tminiprot\tCDS\t10\t100\t.\t+\t0\tParent=MP1\n"
            "contig_1\tmitfi\ttRNA\t200\t220\t.\t-\t.\t"
            "ID=t1;Name=trnA\n"
        ))

    def test_two_panels_when_plastid_substitution_applies(self):
        gff = self._gff()
        metadata = _write(self.tmp, "bin_metadata.json", json.dumps({
            'plastid_canonicalisation': {
                'substitution_applied': True, 'path1_len': 500,
                'lsc_len': 300, 'ir_len': 50, 'ssc_len': 100,
            },
        }))
        panels = om.build_panels(gff, None, metadata)
        labels = [label for label, _, _ in panels]
        self.assertEqual(labels, ['path1', 'path2'])

    def test_single_panel_without_plastid_metadata(self):
        gff = self._gff()
        metadata = _write(self.tmp, "bin_metadata.json", json.dumps({
            'contigs_selected': ['contig_1'],
            'contigs': [{'contig_id': 'contig_1', 'length_bp': 500}],
        }))
        panels = om.build_panels(gff, None, metadata)
        self.assertEqual([label for label, _, _ in panels], ['path1'])

    def test_run_writes_empty_file_when_length_unknown(self):
        gff = _write(self.tmp, "empty.gff", "##gff-version 3\n")
        out = self.tmp / "out.svg"
        om.run(gff, None, None, out)
        self.assertEqual(out.read_text(), '')

    def test_run_writes_svg_with_expected_content(self):
        gff = self._gff()
        out = self.tmp / "out.svg"
        om.run(gff, None, None, out)
        content = out.read_text()
        self.assertIn('<svg', content)
        self.assertIn('data-gene="cox1"', content)
        self.assertIn('data-gene="trnA"', content)

    def test_run_labels_panel_genes_when_panel_given(self):
        gff = self._gff()
        panel = _write(
            self.tmp, "panel.json", json.dumps({"animal_mt": ["cox1"]}))
        out = self.tmp / "out.svg"
        om.run(gff, None, None, out, panel, "animal_mt")
        content = out.read_text()
        self.assertIn('data-barcode="1"', content)
        self.assertIn('data-barcode-label="cox1"', content)
        # The tRNA feature (trnA) is not in the panel and must not be
        # labelled.
        self.assertNotIn('data-barcode-label="trnA"', content)

    def test_unreadable_panel_gives_no_labels_and_still_renders(self):
        gff = self._gff()
        out = self.tmp / "out.svg"
        bad_panel = _write(self.tmp, "bad.json", "{not json")
        om.run(gff, None, None, out, bad_panel, "animal_mt")
        content = out.read_text()
        self.assertIn('<svg', content)
        self.assertNotIn('data-barcode="1"', content)

    def test_target_with_no_panel_entry_gives_no_labels(self):
        gff = self._gff()
        out = self.tmp / "out.svg"
        panel = _write(
            self.tmp, "panel.json", json.dumps({"plant_pt": ["rbcL"]}))
        om.run(gff, None, None, out, panel, "animal_mt")
        content = out.read_text()
        self.assertIn('<svg', content)
        self.assertNotIn('data-barcode="1"', content)


class RenderSvgGeometryTests(unittest.TestCase):

    def _locus(self, start, end, strand='+', kind='CDS', gene='g',
               alt_genes=None, pident=None):
        return {
            'seqid': 'c1', 'source': 'miniprot', 'kind': kind,
            'strand': strand, 'start': start, 'end': end, 'gene': gene,
            'pident': pident, 'qcovhsp': 80.0 if pident else None,
            'bitscore': 100.0 if pident else None,
            'alt_genes': alt_genes or [],
        }

    def test_large_arc_flag_for_wide_feature(self):
        # A feature spanning more than half the genome exercises the
        # large-arc branch of _wedge_path.
        loci = [self._locus(1, 9000)]
        svg = om.render_svg([('path1', loci, 10000)])
        self.assertIn('<svg', svg)

    def test_minimum_wedge_width_applied_to_tiny_feature(self):
        loci = [self._locus(1, 2)]
        svg = om.render_svg([('path1', loci, 1000000)])
        self.assertIn('<path', svg)

    def test_minus_strand_uses_inner_ring(self):
        r_inner, r_outer = om._ring_radii('-')
        self.assertLess(r_outer, om.RING_RADIUS)

    def test_plus_strand_uses_outer_ring(self):
        r_inner, r_outer = om._ring_radii('+')
        self.assertGreater(r_inner, om.RING_RADIUS)

    def test_tooltip_includes_alt_genes_and_scores(self):
        locus = self._locus(
            1, 10, alt_genes=['ALT1'], pident=90.0)
        tooltip = om._tooltip_text(locus)
        self.assertIn('also called: ALT1', tooltip)
        self.assertIn('pident=90.0', tooltip)

    def test_tooltip_without_alt_genes_or_scores(self):
        locus = self._locus(1, 10)
        tooltip = om._tooltip_text(locus)
        self.assertNotIn('also called', tooltip)
        self.assertNotIn('pident', tooltip)

    def test_legend_lists_all_kinds(self):
        legend = om._render_legend(0, 0)
        for kind in om.KIND_COLOURS:
            self.assertIn(kind, legend)

    def test_unknown_kind_uses_default_colour(self):
        loci = [self._locus(1, 10, kind='mystery')]
        svg = om.render_svg([('path1', loci, 1000)])
        self.assertIn(om.DEFAULT_COLOUR, svg)


class RenderBarcodeLabelsTests(unittest.TestCase):

    def _locus(self, start, end, gene='COX1', is_barcode=True):
        return {
            'seqid': 'c1', 'source': 'miniprot', 'kind': 'CDS',
            'strand': '+', 'start': start, 'end': end, 'gene': gene,
            'pident': None, 'qcovhsp': None, 'bitscore': None,
            'alt_genes': [], 'is_barcode': is_barcode,
        }

    def test_no_labels_when_nothing_is_barcode(self):
        loci = [self._locus(1, 10, is_barcode=False)]
        self.assertEqual(om._render_barcode_labels(0, 0, loci, 1000), '')

    def test_single_label_rendered(self):
        loci = [self._locus(1, 10)]
        svg = om._render_barcode_labels(0, 0, loci, 1000)
        self.assertIn('data-barcode-label="COX1"', svg)
        self.assertIn('<line', svg)

    def test_adjacent_labels_nudged_apart(self):
        # Two panel loci close together in angle must not collide:
        # the second label's radius is pushed outward.
        loci = [
            self._locus(1, 10, gene='A'),
            self._locus(11, 20, gene='B'),
        ]
        svg = om._render_barcode_labels(0, 0, loci, 100000)
        texts = re.findall(r'<text x="([\d.-]+)"', svg)
        self.assertEqual(len(texts), 2)
        self.assertNotEqual(texts[0], texts[1])

    def test_anchor_side_follows_angle(self):
        # A locus on the right half of the circle anchors text
        # start-ward; one on the left anchors end-ward, so labels
        # never run back over the ring.
        right = self._locus(1, 10, gene='R')  # near angle -pi/2
        svg_right = om._render_barcode_labels(0, 0, [right], 4)
        self.assertIn('text-anchor="start"', svg_right)


class LegendTextContainmentTests(unittest.TestCase):
    """Task 54 §5.3/§5.4 regression: every legend and label text
    position must lie within the rendered viewBox."""

    def _text_positions(self, svg):
        positions = []
        for m in re.finditer(
                r'<text x="([\d.-]+)" y="([\d.-]+)"[^>]*'
                r'font-size="(\d+)"', svg):
            positions.append(
                (float(m.group(1)), float(m.group(2)), float(m.group(3))))
        return positions

    def test_single_panel_no_labels_fits_viewbox(self):
        loci = [{
            'seqid': 'c1', 'source': 'mitfi', 'kind': 'tRNA',
            'strand': '+', 'start': 1, 'end': 10, 'gene': 'trnA',
            'pident': None, 'qcovhsp': None, 'bitscore': None,
            'alt_genes': [], 'is_barcode': False,
        }]
        svg = om.render_svg([('path1', loci, 1000)])
        vb = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
        width, height = float(vb.group(1)), float(vb.group(2))
        for x, y, fs in self._text_positions(svg):
            self.assertLessEqual(y + fs, height)
            self.assertGreaterEqual(y, 0)
            self.assertLessEqual(x, width)
            self.assertGreaterEqual(x, 0)

    def test_panel_with_several_barcode_labels_fits_viewbox(self):
        loci = []
        length = 20000
        for i, gene in enumerate(
                ['COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ATP6']):
            start = 1 + i * 3000
            loci.append({
                'seqid': 'c1', 'source': 'miniprot', 'kind': 'CDS',
                'strand': '+' if i % 2 == 0 else '-',
                'start': start, 'end': start + 500, 'gene': gene,
                'pident': 90.0, 'qcovhsp': 90.0, 'bitscore': 100.0,
                'alt_genes': [], 'is_barcode': True,
            })
        svg = om.render_svg([('path1', loci, length)])
        vb = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
        width, height = float(vb.group(1)), float(vb.group(2))
        for x, y, fs in self._text_positions(svg):
            self.assertLessEqual(y + fs, height)
            self.assertGreaterEqual(y, 0)
            self.assertLessEqual(x, width)
            self.assertGreaterEqual(x, 0)


class MainTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_main_success(self):
        gff = _write(self.tmp, "a.gff", (
            "contig_1\tminiprot\tmRNA\t10\t100\t.\t+\t.\t"
            "ID=MP1;Target=NC_1_cox1 1 30\n"
            "contig_1\tminiprot\tCDS\t10\t100\t.\t+\t0\tParent=MP1\n"
        ))
        out = self.tmp / "out.svg"
        argv = [
            "organelle_map.py", "--gff", str(gff), "--sample-id", "S1",
            "--out", str(out),
        ]
        old_argv = sys.argv
        sys.argv = argv
        try:
            rc = om.main()
        finally:
            sys.argv = old_argv
        self.assertEqual(rc, 0)
        self.assertIn('<svg', out.read_text())

    def test_main_with_locus_panel_args(self):
        gff = _write(self.tmp, "a.gff", (
            "contig_1\tminiprot\tmRNA\t10\t100\t.\t+\t.\t"
            "ID=MP1;Target=NC_1_cox1 1 30\n"
            "contig_1\tminiprot\tCDS\t10\t100\t.\t+\t0\tParent=MP1\n"
        ))
        panel = _write(
            self.tmp, "panel.json", json.dumps({"animal_mt": ["cox1"]}))
        out = self.tmp / "out.svg"
        argv = [
            "organelle_map.py", "--gff", str(gff), "--sample-id", "S1",
            "--locus-panel", str(panel), "--assembly-target", "animal_mt",
            "--out", str(out),
        ]
        old_argv = sys.argv
        sys.argv = argv
        try:
            rc = om.main()
        finally:
            sys.argv = old_argv
        self.assertEqual(rc, 0)
        self.assertIn('data-barcode="1"', out.read_text())

    def test_main_catches_exceptions_and_exits_0(self):
        out = self.tmp / "out.svg"
        argv = [
            "organelle_map.py", "--gff", str(self.tmp / "missing.gff"),
            "--sample-id", "S1", "--out", str(out),
        ]
        old_argv = sys.argv
        sys.argv = argv
        try:
            rc = om.main()
        finally:
            sys.argv = old_argv
        self.assertEqual(rc, 0)
        self.assertEqual(out.read_text(), '')


if __name__ == "__main__":
    unittest.main()
