#!/usr/bin/env bash
#
# Run the wf5 organelle assembly pipeline on Azure Batch.
#
# Reference data is staged per task from az://refdata-wf5/<version>/.
# Work dir and cache live under az://workdata/wf5/ and az://cache/wf5/.
#
# Usage:
#   deploy/azure/run-wf5-azure.sh --stub --refs-version v2026.09_1
#   deploy/azure/run-wf5-azure.sh --integration --refs-version v2026.09_1
#   deploy/azure/run-wf5-azure.sh --refs-version v2026.09_1 \
#       --samplesheet samples.csv [--data-dir /path/to/reads] \
#       [--outdir out] [--pool wf5] [--polish] [-resume] [--trace] [--yes] \
#       [-- extra nextflow args]
#
# Omit --data-dir when the samplesheet lists az:// read URIs.

set -euo pipefail

NXF_VERSION="25.04.0"
ENV_FILE=".env.azure"
REFDATA_CONTAINER="refdata-wf5"
WORK_CONTAINER="workdata"
WORK_PREFIX="wf5"
CACHE_CONTAINER="cache"
CACHE_PREFIX="wf5"
FIXTURES_CONTAINER="test-data"
DEFAULT_POOL="wf5_test"
PROD_POOL="wf5"
DEFAULT_FIXTURES_VERSION="v2026.07"
INTEGRATION_SAMPLESHEET="tests/integration/samples.csv"
INTEGRATION_ASSERTIONS="tests/integration/assertions.sh"
REQUIRED_VARS=(
    AZURE_BATCH_ACCOUNT_NAME
    AZURE_BATCH_ENDPOINT
    AZURE_BATCH_ACCESS_KEY
    STORAGE_ACCOUNT_STD
    AZURE_STORAGE_ACCOUNT_KEY
)

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

SAMPLESHEET=""
DATA_DIR=""
OUTDIR=""
REFS_VERSION=""
FIXTURES_VERSION="${FIXTURES_VERSION:-${DEFAULT_FIXTURES_VERSION}}"
POOL="${DEFAULT_POOL}"
POLISH=""
RESUME=""
TRACE=false
YES=false
EXTRA_ARGS=()
MODE="custom"
RUN_ID="$(date +"%Y%m%d_%H%M%S")_$$"

die() {
    echo -e "${RED}ERROR: $*${NC}" >&2
    exit 1
}

usage() {
    sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-1}"
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
            --data-dir) DATA_DIR="$2"; shift 2 ;;
            --outdir) OUTDIR="$2"; shift 2 ;;
            --refs-version) REFS_VERSION="$2"; shift 2 ;;
            --fixtures-version) FIXTURES_VERSION="$2"; shift 2 ;;
            --pool) POOL="$2"; shift 2 ;;
            --polish) POLISH="--polish"; shift ;;
            -resume|--resume) RESUME="-resume"; shift ;;
            --trace) TRACE=true; shift ;;
            --stub) MODE="stub"; shift ;;
            --integration) MODE="integration"; shift ;;
            --yes|-y) YES=true; shift ;;
            --) shift; EXTRA_ARGS=("$@"); break ;;
            -h|--help) usage 0 ;;
            *) echo "Unknown argument: $1" >&2; usage 1 ;;
        esac
    done
}

load_env() {
    [[ -f "${ENV_FILE}" ]] \
        || die "${ENV_FILE} not found (run from the repository root)"
    set -a
    # shellcheck disable=SC1090
    source "${ENV_FILE}"
    set +a

    local var
    for var in "${REQUIRED_VARS[@]}"; do
        [[ -n "${!var:-}" ]] || die "${var} is not set in ${ENV_FILE}"
    done
}

validate_args() {
    [[ -n "${REFS_VERSION}" ]] || die "--refs-version is required"

    case "${POOL}" in
        "${DEFAULT_POOL}") CACHE_ENV="test" ;;
        "${PROD_POOL}") CACHE_ENV="prod" ;;
        *) die "--pool must be '${DEFAULT_POOL}' or '${PROD_POOL}', got '${POOL}'" ;;
    esac

    case "${MODE}" in
        stub)
            [[ -z "${SAMPLESHEET}" && -z "${DATA_DIR}" ]] \
                || die "--stub uses the built-in stub fixtures; omit --samplesheet/--data-dir"
            ;;
        integration)
            [[ -z "${SAMPLESHEET}" && -z "${DATA_DIR}" ]] \
                || die "--integration uses the built-in fixtures; omit --samplesheet/--data-dir"
            SAMPLESHEET="${INTEGRATION_SAMPLESHEET}"
            DATA_DIR="az://${FIXTURES_CONTAINER}/wf5/${FIXTURES_VERSION}"
            ;;
        custom)
            [[ -n "${SAMPLESHEET}" ]] || die "--samplesheet is required"
            [[ -f "${SAMPLESHEET}" ]] || die "samplesheet not found: ${SAMPLESHEET}"
            if [[ -n "${DATA_DIR}" && "${DATA_DIR}" != az://* ]]; then
                [[ -d "${DATA_DIR}" ]] || die "--data-dir not found: ${DATA_DIR}"
            fi
            ;;
    esac

    if [[ -z "${OUTDIR}" ]]; then
        OUTDIR="output/${RUN_ID}"
        [[ "${MODE}" == "custom" ]] || OUTDIR="output/${MODE}_${RUN_ID}"
    fi
}

batch_az() {
    az batch "$@" \
        --account-name "${AZURE_BATCH_ACCOUNT_NAME}" \
        --account-endpoint "${AZURE_BATCH_ENDPOINT#https://}" \
        --account-key "${AZURE_BATCH_ACCESS_KEY}"
}

preflight_pool() {
    local state
    state=$(batch_az pool show --pool-id "${POOL}" --query state -o tsv 2>/dev/null) \
        || die "Batch pool '${POOL}' does not exist or is unreachable"
    [[ "${state}" == "active" ]] \
        || die "Batch pool '${POOL}' is in state '${state}', not 'active'"
}

preflight_storage() {
    local exists
    exists=$(az storage container exists \
        --account-name "${STORAGE_ACCOUNT_STD}" \
        --account-key "${AZURE_STORAGE_ACCOUNT_KEY}" \
        --name "${WORK_CONTAINER}" --query exists -o tsv 2>/dev/null) \
        || die "cannot reach storage account ${STORAGE_ACCOUNT_STD} with the configured key"
    [[ "${exists}" == "true" ]] \
        || die "container '${WORK_CONTAINER}' not found in ${STORAGE_ACCOUNT_STD}"
}

# Mode A uploads the listed reads over this machine's link, so show the
# cost up front. Reads are the sheet's 'reads' column, relative to data dir.
print_upload_size() {
    [[ -n "${DATA_DIR}" && "${DATA_DIR}" != az://* ]] || return 0
    local total
    total=$(python3 - "${SAMPLESHEET}" "${DATA_DIR}" <<'PYEOF'
import csv
import sys
from pathlib import Path

sheet, data_dir = sys.argv[1], Path(sys.argv[2])
total = 0
for row in csv.DictReader(open(sheet)):
    for name in row["reads"].replace(";", " ").replace(",", " ").split():
        path = Path(name)
        if not path.is_absolute():
            path = data_dir / path
        total += path.stat().st_size if path.exists() else 0
print(f"{total / 1e9:.2f} GB")
PYEOF
)
    echo "Reads to upload: ${total} (from ${DATA_DIR})"
}

confirm() {
    [[ "${YES}" == true ]] && return 0
    local reply
    read -r -p "Continue with workflow execution? (yes/no): " reply
    [[ "${reply}" == "yes" ]] || { echo "Execution cancelled"; exit 0; }
}

build_command() {
    local az_refs="az://${REFDATA_CONTAINER}/${REFS_VERSION}"
    local profile="azure"
    [[ "${MODE}" != "stub" ]] || profile="azure,stub"
    NF_ARGS=(
        run main.nf
        -profile "${profile}"
        --outdir "${OUTDIR}"
        --organelle_refs "${az_refs}/recruit"
        --protein_panel "${az_refs}/proteins"
        --blast_db "${az_refs}/validate"
        --annotate_refs "${az_refs}/annotate"
        --refs_manifest "${az_refs}/manifest.json"
    )
    [[ -z "${POLISH}" ]] || NF_ARGS+=("${POLISH}")
    [[ -z "${RESUME}" ]] || NF_ARGS+=("${RESUME}")
    # Trace collection needs `ps` inside each task container
    [[ "${TRACE}" == false ]] || NF_ARGS+=(-with-trace "${OUTDIR}/trace.txt")

    case "${MODE}" in
        stub)
            NF_ARGS+=(
                -stub-run
                --samplesheet "tests/stub/samplesheet.csv"
                --data_dir "tests/stub/data"
            )
            ;;
        integration)
            NF_ARGS+=(
                --samplesheet "${SAMPLESHEET}"
                --data_dir "${DATA_DIR}"
                --publish_intermediates true
            )
            ;;
        custom)
            NF_ARGS+=(--samplesheet "${SAMPLESHEET}")
            [[ -z "${DATA_DIR}" ]] || NF_ARGS+=(--data_dir "${DATA_DIR}")
            ;;
    esac
    NF_ARGS+=("${EXTRA_ARGS[@]}")
}

main() {
    cd "$(dirname "${BASH_SOURCE[0]}")/../.."
    parse_args "$@"
    load_env
    validate_args
    preflight_pool
    preflight_storage

    export NXF_VER="${NXF_VERSION}"
    export WF5_AZURE_POOL="${POOL}"
    export WF5_AZURE_WORKDIR="az://${WORK_CONTAINER}/${WORK_PREFIX}/${RUN_ID}"
    export WF5_AZURE_CACHE="az://${CACHE_CONTAINER}/${CACHE_PREFIX}/${CACHE_ENV}"
    build_command

    echo -e "${YELLOW}=== wf5 Azure Batch run ===${NC}"
    echo "Mode:            ${MODE}"
    echo "Pool:            ${POOL}"
    echo "Refs version:    ${REFS_VERSION} (az://${REFDATA_CONTAINER}/${REFS_VERSION}/)"
    echo "Samplesheet:     ${SAMPLESHEET:-(stub fixtures)}"
    echo "Data dir:        ${DATA_DIR:-(none; sheet lists az:// URIs)}"
    echo "Outdir:          ${OUTDIR}"
    echo "Work dir:        ${WF5_AZURE_WORKDIR}"
    echo "Cache:           ${WF5_AZURE_CACHE}"
    echo "Nextflow:        ${NXF_VER}"
    echo "Polish:          ${POLISH:-off}"
    echo "Resume:          ${RESUME:-false}"
    print_upload_size
    echo ""
    confirm

    mkdir -p "${OUTDIR}"
    local exit_code=0
    nextflow "${NF_ARGS[@]}" || exit_code=$?

    if [[ ${exit_code} -ne 0 ]]; then
        echo -e "${RED}=== Workflow failed (exit ${exit_code}) ===${NC}"
        echo "Pull task logs: source deploy/azure/batch-helpers.sh && az_fetch_nf_logs .nextflow.log ./logs"
        exit "${exit_code}"
    fi

    echo -e "${GREEN}=== Workflow completed ===${NC}"
    echo "Work dir: ${WF5_AZURE_WORKDIR}"
    echo "Run report: ${OUTDIR}/run-report.html"

    if [[ "${MODE}" == "integration" ]]; then
        NF_LOG_FILE=.nextflow.log bash "${INTEGRATION_ASSERTIONS}" "${OUTDIR}"
    fi
}

main "$@"
