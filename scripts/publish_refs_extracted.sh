#!/usr/bin/env bash
# Publishes an unpacked reference bundle to Azure blob, beside the tarball:
#   refdata-wf5/<version>/{recruit,validate,proteins,annotate,manifest.json}
# Azure Batch runs stage these directories straight from blob.
#
# Bundles are immutable: refuses to run if the extracted form already
# exists for <version>. manifest.json is uploaded last and acts as the
# commit marker. The upload is then checked against the manifest's
# artefacts list, and fails loudly if anything is missing.
#
# Needs .env.azure (STORAGE_ACCOUNT_STD, AZURE_STORAGE_ACCOUNT_KEY).
#
# Usage:
#   bash scripts/publish_refs_extracted.sh v2026.10_1

set -euo pipefail

VERSION="${1:?Usage: publish_refs_extracted.sh <version>  e.g. v2026.10_1}"
REFDATA_CONTAINER="${REFDATA_CONTAINER:-refdata-wf5}"
BUNDLE_DIR="refs/${VERSION}"
ENV_FILE=".env.azure"

if [[ ! -f "${BUNDLE_DIR}/manifest.json" ]]; then
    echo "ERROR: ${BUNDLE_DIR}/manifest.json not found" \
        "(run scripts/fetch_refs.sh ${VERSION} first)" >&2
    exit 1
fi

set -a
# shellcheck disable=SC1090
source "${ENV_FILE}"
set +a
: "${STORAGE_ACCOUNT_STD:?not set in ${ENV_FILE}}"
: "${AZURE_STORAGE_ACCOUNT_KEY:?not set in ${ENV_FILE}}"

AZ_AUTH=(
    --account-name "${STORAGE_ACCOUNT_STD}"
    --account-key "${AZURE_STORAGE_ACCOUNT_KEY}"
)

list_blobs() {
    az storage blob list "${AZ_AUTH[@]}" \
        --container-name "${REFDATA_CONTAINER}" \
        --prefix "$1" --query '[].name' -o tsv
}

EXISTING=$(list_blobs "${VERSION}/" | grep -v "^${VERSION}/refs.tar.gz$" || true)
if [[ -n "${EXISTING}" ]]; then
    echo "ERROR: extracted bundle already exists for ${VERSION}; bundles" \
        "are immutable. Publish under a new version." >&2
    echo "${EXISTING}" | head -5 >&2
    exit 1
fi

COMPONENTS=$(python3 - "${BUNDLE_DIR}/manifest.json" <<'PYEOF'
import json
import sys

manifest = json.load(open(sys.argv[1]))
print("\n".join(sorted({p.split("/", 1)[0] for p in manifest["artefacts"]})))
PYEOF
)

for component in ${COMPONENTS}; do
    echo "Uploading ${component}/ ..."
    az storage blob upload-batch "${AZ_AUTH[@]}" \
        --destination "${REFDATA_CONTAINER}" \
        --destination-path "${VERSION}/${component}" \
        --source "${BUNDLE_DIR}/${component}" \
        --overwrite false --only-show-errors > /dev/null
done

echo "Uploading manifest.json ..."
az storage blob upload "${AZ_AUTH[@]}" \
    --container-name "${REFDATA_CONTAINER}" \
    --name "${VERSION}/manifest.json" \
    --file "${BUNDLE_DIR}/manifest.json" \
    --overwrite false --only-show-errors > /dev/null

echo "Verifying blob listing against manifest.json artefacts..."
list_blobs "${VERSION}/" | sort > "${TMPDIR:-/tmp}/refs_blobs.$$"
MISSING=$(python3 - "${BUNDLE_DIR}/manifest.json" "${VERSION}" \
    "${TMPDIR:-/tmp}/refs_blobs.$$" <<'PYEOF'
import json
import sys

manifest = json.load(open(sys.argv[1]))
version = sys.argv[2]
blobs = set(open(sys.argv[3]).read().split("\n"))
expected = {f"{version}/{p}" for p in manifest["artefacts"]}
expected.add(f"{version}/manifest.json")
print("\n".join(sorted(expected - blobs)))
PYEOF
)
rm -f "${TMPDIR:-/tmp}/refs_blobs.$$"

if [[ -n "${MISSING}" ]]; then
    echo "ERROR: published bundle is incomplete; missing:" >&2
    echo "${MISSING}" | head -20 >&2
    exit 1
fi

echo "Published az://${REFDATA_CONTAINER}/${VERSION}/ (complete)"
