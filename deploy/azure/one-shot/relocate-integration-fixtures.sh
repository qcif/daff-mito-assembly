#!/usr/bin/env bash
# One-shot: move integration-fixtures/<version>/* → integration-fixtures/wf5/<version>/*
# on daffstandard, to match the workflow-scoped layout in spec/04a-azure-blob-storage.md.
#
# Prereqs:
#   - az login as an account with Storage Blob Data Contributor on daffstandard.
#
# Usage:
#   bash deploy/azure/relocate-integration-fixtures.sh [version]
#
# Defaults to v2026.07. Verifies SHA256 for each copied blob against
# tests/integration/fetched.sha256 before deleting the originals.

set -euo pipefail

source .env.azure

ACCOUNT=$STORAGE_ACCOUNT_STD
CONTAINER=$STORAGE_CONTAINER_FIXTURES
VERSION="${1:-v2026.07}"
SHA_FILE="tests/integration/fetched.sha256"

AZ=(--account-name "$ACCOUNT" --account-key "$AZURE_STORAGE_ACCOUNT_KEY")

echo "== 1. List blobs under ${CONTAINER}/${VERSION}/ =="
mapfile -t BLOBS < <(
    az storage blob list "${AZ[@]}" \
        --container-name "$CONTAINER" \
        --prefix "${VERSION}/" \
        --query '[].name' -o tsv
)

if [[ ${#BLOBS[@]} -eq 0 ]]; then
    echo "No blobs under ${VERSION}/ — nothing to relocate. Already migrated?"
    exit 0
fi

printf '  %s\n' "${BLOBS[@]}"

echo
echo "== 2. Server-side copy to wf5/${VERSION}/ =="
for src in "${BLOBS[@]}"; do
    # Strip leading "<version>/" prefix, prepend "wf5/<version>/".
    rel="${src#${VERSION}/}"
    dst="wf5/${VERSION}/${rel}"
    echo "  ${src}  →  ${dst}"
    az storage blob copy start "${AZ[@]}" \
        --source-container "$CONTAINER" \
        --source-blob "$src" \
        --destination-container "$CONTAINER" \
        --destination-blob "$dst" >/dev/null
done

echo
echo "== 3. Wait for server-side copies to complete =="
for src in "${BLOBS[@]}"; do
    rel="${src#${VERSION}/}"
    dst="wf5/${VERSION}/${rel}"
    while :; do
        status=$(
            az storage blob show "${AZ[@]}" \
                --container-name "$CONTAINER" \
                --name "$dst" \
                --query 'properties.copy.status' -o tsv
        )
        case "$status" in
            success) echo "  OK   $dst"; break ;;
            pending) sleep 2 ;;
            *)       echo "  FAIL $dst (status=$status)"; exit 1 ;;
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
        rel="${src#${VERSION}/}"
        dst="wf5/${VERSION}/${rel}"
        local_path="${TMPDIR}/${rel}"
        mkdir -p "$(dirname "$local_path")"
        az storage blob download "${AZ[@]}" \
            --container-name "$CONTAINER" \
            --name "$dst" \
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
echo "== 5. Delete originals under ${VERSION}/ =="
for src in "${BLOBS[@]}"; do
    echo "  rm $src"
    az storage blob delete "${AZ[@]}" \
        --container-name "$CONTAINER" \
        --name "$src" >/dev/null
done

echo
echo "Done. New layout:"
az storage blob list "${AZ[@]}" \
    --container-name "$CONTAINER" \
    --prefix "wf5/${VERSION}/" \
    --query '[].name' -o tsv
