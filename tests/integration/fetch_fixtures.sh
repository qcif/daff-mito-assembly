#!/usr/bin/env bash
# Downloads Tier 2 integration fixtures from Azure blob storage (public access).
#
# Usage:
#   bash tests/integration/fetch_fixtures.sh

set -euo pipefail

STORAGE_ACCOUNT="${AZURE_STORAGE_ACCOUNT:-daffstandard}"
FIXTURES_CONTAINER="${FIXTURES_CONTAINER:-test-data}"
FIXTURES_VERSION="${FIXTURES_VERSION:-v2026.10}"
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

# Deterministic corrupt fixture: a truncated gzip stream, which the intake
# check must reject (INT-CORRUPT-01). Generated after checksum verification
# since it is not a pinned fixture.
CORRUPT_BYTES=1000000
head -c "$CORRUPT_BYTES" "${FIXTURES_DIR}/animal_mt.fastq.gz" \
    > "${FIXTURES_DIR}/INT-CORRUPT-01.fastq.gz"

# Deterministic duplicate-ID fixture: the first N records of animal_mt
# followed by the same N records again (INT-DUPID-01). N is small enough
# that the sample stops at the coverage gate.
DUPID_RECORDS=10
zcat "${FIXTURES_DIR}/animal_mt.fastq.gz" \
    | head -n $(( DUPID_RECORDS * 4 )) \
    > "${FIXTURES_DIR}/dupid_head.fastq"
cat "${FIXTURES_DIR}/dupid_head.fastq" "${FIXTURES_DIR}/dupid_head.fastq" \
    | gzip -n > "${FIXTURES_DIR}/INT-DUPID-01.fastq.gz"
rm "${FIXTURES_DIR}/dupid_head.fastq"

echo "Fixtures ready in $FIXTURES_DIR"
