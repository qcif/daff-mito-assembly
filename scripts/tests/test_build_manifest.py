"""Unit tests for scripts/refdata/build_manifest.py (task 64 §3.1):
the protein panel's `label_check` provenance is folded into the
manifest."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[2]
    / "scripts" / "refdata" / "build_manifest.py"
)
_spec = importlib.util.spec_from_file_location("build_manifest", MODULE_PATH)
bm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bm)

LABEL_CHECK = {"plant_pt": {
    "method": "m", "iterations": 2, "no_hits": 3, "dropped_genes": [],
    "excluded": [{
        "accession": "NC_034951.1", "labelled_gene": "ycf3",
        "best_other_gene": "clpP",
        "own_bitscore": 41.2, "other_bitscore": 389.0}]}}


class TestManifestLabelCheck(unittest.TestCase):

    def test_label_check_folded_into_protein_panel(self):
        with tempfile.TemporaryDirectory() as tmp:
            bundle = Path(tmp) / "v9"
            (bundle / "proteins").mkdir(parents=True)
            (bundle / "proteins" / "provenance.json").write_text(
                json.dumps({"refseq_release": "237",
                            "label_check": LABEL_CHECK}))
            manifest = bm.build_manifest(bundle, "db")
        self.assertEqual(
            manifest["protein_panel"]["label_check"], LABEL_CHECK)
        self.assertNotIn(
            "proteins/provenance.json", manifest["artefacts"])


if __name__ == "__main__":
    unittest.main()
