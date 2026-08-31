#!/usr/bin/env bash
# Downloads and unpacks the reference bundle from Azure blob (public access).
#
# Usage:
#   bash scripts/fetch_refs.sh v2026.07

set -euo pipefail

VERSION="${1:?Usage: fetch_refs.sh <version>  e.g. v2026.07}"
STORAGE_ACCOUNT="${AZURE_STORAGE_ACCOUNT:-daffstandard}"
REFDATA_CONTAINER="${REFDATA_CONTAINER:-refdata-wf5}"
BASE_URL="https://${STORAGE_ACCOUNT}.blob.core.windows.net/${REFDATA_CONTAINER}"
TARBALL="refs-${VERSION}.tar.gz"

if [ -d "refs/${VERSION}" ]; then
    echo "refs/${VERSION} already present — skipping download."
    exit 0
fi

mkdir -p refs

echo "Downloading ${VERSION} refdata bundle..."
curl -fsSL "${BASE_URL}/${VERSION}/refs.tar.gz" -o "${TARBALL}"

echo "Verifying SHA256..."
sha256sum -c <(sed "s|refs.tar.gz|${TARBALL}|" "scripts/refs-${VERSION}.sha256")

echo "Unpacking..."
tar -xzf "${TARBALL}" -C refs
rm "${TARBALL}"

echo "Verifying bundle completeness against manifest.json..."
MISSING=$(python3 - "${VERSION}" <<'PYEOF'
import json
import sys
from pathlib import Path

version = sys.argv[1]
bundle_root = Path("refs") / version
manifest = json.loads((bundle_root / "manifest.json").read_text())
components = sorted({p.split("/", 1)[0] for p in manifest["artefacts"]})
missing = [c for c in components if not (bundle_root / c).is_dir()]
print("\n".join(missing))
PYEOF
)

if [ -n "${MISSING}" ]; then
    echo "ERROR: refs/${VERSION} is missing components listed in" \
        "manifest.json: $(echo "${MISSING}" | tr '\n' ' ')" >&2
    exit 1
fi

echo "Refdata ready in refs/${VERSION}/"
