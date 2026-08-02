#!/usr/bin/env bash
# One-shot: move integration-fixtures/wf5/<version>/* → test-data/wf5/<version>/*
# on daffstandard, consolidating fixtures into the shared test-data container.
#
# Prereqs:
#   - .env.azure with AZURE_STORAGE_ACCOUNT_KEY for daffstandard.
#   - Destination container must already exist with the same public access
#     level as the source (blob), since fetch_fixtures.sh reads anonymously.
#
# Usage:
#   bash deploy/azure/relocate-fixtures-container.sh [version]
#
# Defaults to v2026.07. Verifies SHA256 for each copied blob against
# tests/integration/fetched.sha256 before deleting the originals.

set -euo pipefail

source .env.azure

ACCOUNT=$STORAGE_ACCOUNT_STD
SRC_CONTAINER=integration-fixtures
DST_CONTAINER=test-data
VERSION="${1:-v2026.07}"
PREFIX="wf5/${VERSION}"
SHA_FILE="tests/integration/fetched.sha256"

AZ=(--account-name "$ACCOUNT" --account-key "$AZURE_STORAGE_ACCOUNT_KEY")

echo "== 1. List blobs under ${SRC_CONTAINER}/${PREFIX}/ =="
mapfile -t BLOBS < <(
    az storage blob list "${AZ[@]}" \
        --container-name "$SRC_CONTAINER" \
        --prefix "${PREFIX}/" \
        --query '[].name' -o tsv
)

if [[ ${#BLOBS[@]} -eq 0 ]]; then
    echo "No blobs under ${PREFIX}/ — nothing to relocate. Already migrated?"
    exit 0
fi

printf '  %s\n' "${BLOBS[@]}"

echo
echo "== 2. Server-side copy to ${DST_CONTAINER}/ (same blob names) =="
for src in "${BLOBS[@]}"; do
    echo "  ${SRC_CONTAINER}/${src}  →  ${DST_CONTAINER}/${src}"
    az storage blob copy start "${AZ[@]}" \
        --source-container "$SRC_CONTAINER" \
        --source-blob "$src" \
        --destination-container "$DST_CONTAINER" \
        --destination-blob "$src" >/dev/null
done

echo
echo "== 3. Wait for server-side copies to complete =="
for src in "${BLOBS[@]}"; do
    while :; do
        status=$(
            az storage blob show "${AZ[@]}" \
                --container-name "$DST_CONTAINER" \
                --name "$src" \
                --query 'properties.copy.status' -o tsv
        )
        case "$status" in
            success) echo "  OK   $src"; break ;;
            pending) sleep 2 ;;
            *)       echo "  FAIL $src (status=$status)"; exit 1 ;;
        esac
    done
done

echo
echo "== 4. Verify SHA256 of copied blobs against ${SHA_FILE} =="
if [[ ! -f "$SHA_FILE" ]]; then
    echo "WARN: $SHA_FILE not found — skipping content verification." >&2
else
    TMPDIR=$(mktemp -d)
    trap 'rm -rf "$TMPDIR"' EXIT
    for src in "${BLOBS[@]}"; do
        rel="${src#${PREFIX}/}"
        local_path="${TMPDIR}/${rel}"
        mkdir -p "$(dirname "$local_path")"
        az storage blob download "${AZ[@]}" \
            --container-name "$DST_CONTAINER" \
            --name "$src" \
            --file "$local_path" >/dev/null
        expected=$(awk -v k="$rel" '$2 == k {print $1}' "$SHA_FILE")
        actual=$(sha256sum "$local_path" | awk '{print $1}')
        if [[ -z "$expected" ]]; then
            echo "  SKIP $rel (not in $SHA_FILE)"
        elif [[ "$expected" == "$actual" ]]; then
            echo "  OK   $rel"
        else
            echo "  FAIL $rel  expected=$expected actual=$actual" >&2
            exit 1
        fi
    done
fi

echo
echo "== 5. Delete originals from ${SRC_CONTAINER} =="
for src in "${BLOBS[@]}"; do
    echo "  rm ${SRC_CONTAINER}/${src}"
    az storage blob delete "${AZ[@]}" \
        --container-name "$SRC_CONTAINER" \
        --name "$src" >/dev/null
done

echo
echo "Done. New layout:"
az storage blob list "${AZ[@]}" \
    --container-name "$DST_CONTAINER" \
    --prefix "${PREFIX}/" \
    --query '[].name' -o tsv
