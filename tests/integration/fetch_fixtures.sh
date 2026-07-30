#!/usr/bin/env bash
# Downloads Tier 2 integration fixtures from Azure blob storage (public access).
#
# Usage:
#   bash tests/integration/fetch_fixtures.sh

set -euo pipefail

STORAGE_ACCOUNT="${AZURE_STORAGE_ACCOUNT:-daffstandard}"
FIXTURES_CONTAINER="${FIXTURES_CONTAINER:-integration-fixtures}"
FIXTURES_VERSION="${FIXTURES_VERSION:-v2026.07}"
FIXTURES_DIR="tests/integration/fetched"
BASE_URL="https://${STORAGE_ACCOUNT}.blob.core.windows.net/${FIXTURES_CONTAINER}/wf5/${FIXTURES_VERSION}"

mkdir -p "$FIXTURES_DIR"

echo "Fetching Tier 2 fixtures from Azure blob storage..."

for f in animal_mt plant_pt plant_mt; do
    echo "  Downloading ${f}.fastq.gz ..."
    curl -fsSL "${BASE_URL}/${f}.fastq.gz" -o "${FIXTURES_DIR}/${f}.fastq.gz"
done

echo "Verifying SHA256 checksums..."
python3 tests/integration/verify_shas.py \
    "$FIXTURES_DIR" \
    tests/integration/fetched.sha256

echo "Fixtures ready in $FIXTURES_DIR"
