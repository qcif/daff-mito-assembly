#!/usr/bin/env python3
"""
Emit refs/<version>/manifest.json — spec §4.4.

`manifest.json` is the single source of truth for reference provenance:
its version string is emitted into per-run metadata so every result is
traceable to the exact bundle used (CONSTITUTION rule 10). This script
walks a built bundle, digests every artefact, and folds in the
provenance sidecars the individual build scripts leave behind:

  * `validate/provenance.json`   — RefSeq release + download date
  * `recruit/plant_mt.masking.json` — the task 28 §3.2 masking
    parameters actually applied, panel-wide masked fraction, the worst
    per-genome offenders, and any genome dropped by the max-masked
    guard.

An auditor reading a result six months from now must be able to tell
that the `plant_mt` panel was masked, by what rule, and how much was
removed (task 28 §3.4).

Usage:
    python3 build_manifest.py --bundle refs/v2026.08 [--getorganelle-db DIR]
"""

import argparse
import hashlib
import json
import sys
from datetime import date
from pathlib import Path

CHUNK = 1 << 20

# Directories walked for artefact digests, in manifest order.
BUNDLE_DIRS = ["recruit", "validate", "proteins"]

# Sidecars folded into the manifest rather than digested as artefacts.
PROVENANCE_SIDECAR = Path("validate/provenance.json")
MASKING_SIDECAR = Path("recruit/plant_mt.masking.json")
PROTEIN_PANEL_SIDECAR = Path("proteins/provenance.json")

GETORGANELLE_URL = "https://github.com/Kinggerm/GetOrganelleDB"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", required=True,
                   help="built bundle root, e.g. refs/v2026.08")
    p.add_argument("--getorganelle-db", default="",
                   help="GetOrganelleDB release tag or checkout path")
    return p.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(CHUNK), b''):
            digest.update(block)
    return digest.hexdigest()


def digest_tree(bundle: Path) -> dict:
    """{relative path: {bytes, sha256}} for every artefact in the bundle."""
    artefacts = {}
    for subdir in BUNDLE_DIRS:
        root = bundle / subdir
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(bundle)
            if rel in (PROVENANCE_SIDECAR, MASKING_SIDECAR,
                       PROTEIN_PANEL_SIDECAR):
                continue
            artefacts[str(rel)] = {
                'bytes': path.stat().st_size,
                'sha256': sha256(path),
            }
    return artefacts


def read_json(path: Path):
    """Parse `path`, or return None if it is absent."""
    if not path.is_file():
        print(f"[build_manifest] WARNING: no {path}", file=sys.stderr)
        return None
    return json.loads(path.read_text())


def recruit_panels(masking: dict, getorganelle_db: str) -> dict:
    """
    Per-target recruitment-panel provenance.

    `plant_pt` and `animal_mt` are GetOrganelle seed copies. `plant_mt`
    is derived from RefSeq and masked, so it carries the full masking
    record (task 28 §3.4).
    """
    panels = {
        target: {
            'source': 'GetOrganelleDB SeedDatabase',
            'source_url': GETORGANELLE_URL,
            'release': getorganelle_db,
            'derivation': 'copied verbatim',
        }
        for target in ('plant_pt', 'animal_mt')
    }

    plant_mt = {
        'source': 'RefSeq mitochondrion, Viridiplantae subset',
        'source_file': str(PROVENANCE_SIDECAR.parent
                           / 'refseq_mt_viridiplantae.fa'),
        'derivation': (
            'plastid-derived (NUPT) regions masked to N against the '
            'plant_pt panel — see scripts/refdata/mask_panel.py'
        ),
    }
    if masking is None:
        plant_mt['masking'] = None
        print(
            "[build_manifest] WARNING: plant_mt masking report missing; "
            "manifest cannot record how the panel was derived",
            file=sys.stderr,
        )
    else:
        plant_mt['masking'] = masking
    panels['plant_mt'] = plant_mt

    return panels


def build_manifest(bundle: Path, getorganelle_db: str) -> dict:
    provenance = read_json(bundle / PROVENANCE_SIDECAR)
    masking = read_json(bundle / MASKING_SIDECAR)
    protein_panel = read_json(bundle / PROTEIN_PANEL_SIDECAR)

    return {
        '$schema': 'wf5/refs-manifest/v1',
        'version': bundle.name,
        'generated_at': date.today().isoformat(),
        'refseq': provenance,
        'recruit_panels': recruit_panels(masking, getorganelle_db),
        'protein_panel': protein_panel,
        'artefacts': digest_tree(bundle),
    }


def main():
    args = parse_args()
    bundle = Path(args.bundle)
    if not bundle.is_dir():
        print(f"ERROR: no such bundle: {bundle}", file=sys.stderr)
        sys.exit(1)

    manifest = build_manifest(bundle, args.getorganelle_db)
    out = bundle / "manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"[build_manifest] {out} — {len(manifest['artefacts'])} artefacts",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
