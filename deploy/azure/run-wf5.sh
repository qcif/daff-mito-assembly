#!/bin/bash
#
# Run the wf5 organelle assembly pipeline on Azure Batch.
#
# Reference data is pre-staged to compute nodes by the start task (setup.sh)
# under /mnt/refdata/v<YYYY.MM>/ per plan.md §4.4.
#
# Usage:
#   ./deploy/azure/run-wf5.sh \
#       --samplesheet /path/to/samples.csv \
#       --data-dir /path/to/reads \
#       --outdir /path/to/output
#

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

PID=$$
RUN_ID="$(date +"%Y%m%d_%H%M%S")_$PID"

# Defaults
SAMPLESHEET="samples.csv"
DATA_DIR=""
OUTDIR="output/$RUN_ID"
RESUME=""

# Reference bundle path on the compute node (staged by setup.sh).
# TODO: bump when a new versioned bundle is built (plan.md §4.4).
REFDATA_VERSION="v2026.06"
REFDATA_ROOT="/mnt/refdata/${REFDATA_VERSION}"

# Per-artefact paths within the staged bundle
KINGDOM_REFS="${REFDATA_ROOT}/recruit"
BLAST_DB="${REFDATA_ROOT}/validate"
LOCUS_PANEL="${REFDATA_ROOT}/proteins"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samplesheet)
            SAMPLESHEET="$2"
            shift 2
            ;;
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --polish)
            POLISH_FLAG="--polish"
            shift
            ;;
        -resume)
            RESUME="-resume"
            shift
            ;;
        *)
            echo -e "${RED}ERROR: Unknown argument: $1${NC}"
            echo "Usage: $0 --samplesheet <csv> --data-dir <dir> [--outdir <dir>] [--polish] [-resume]"
            exit 1
            ;;
    esac
done

if [[ -z "$DATA_DIR" ]]; then
    echo -e "${RED}ERROR: --data-dir is required${NC}"
    exit 1
fi

if [[ ! -f .env.azure ]]; then
    echo -e "${RED}ERROR: .env.azure not found${NC}"
    echo "Please run this script from the repository root directory"
    exit 1
fi

echo -e "${YELLOW}Loading Azure credentials from .env.azure${NC}"
set -a
source .env.azure
set +a

if [[ -z "${AZURE_STORAGE_ACCOUNT_KEY:-}" ]]; then
    echo -e "${RED}ERROR: AZURE_STORAGE_ACCOUNT_KEY not set in .env.azure${NC}"
    exit 1
fi

if [[ -z "${AZURE_BATCH_ACCESS_KEY:-}" ]]; then
    echo -e "${RED}ERROR: AZURE_BATCH_ACCESS_KEY not set in .env.azure${NC}"
    exit 1
fi

echo ""
echo -e "${YELLOW}=== wf5 Azure Batch Configuration ===${NC}"
echo "Samplesheet:     $SAMPLESHEET"
echo "Data dir:        $DATA_DIR"
echo "Output dir:      $OUTDIR"
echo "Reference bundle: $REFDATA_ROOT (staged on Azure Batch nodes)"
echo "  Kingdom refs:  $KINGDOM_REFS"
echo "  BLAST DB:      $BLAST_DB"
echo "  Locus panel:   $LOCUS_PANEL"
echo "Profile:         azure"
echo "Polish:          ${POLISH_FLAG:-off}"
echo "Resume:          ${RESUME:-false}"
echo ""

read -p "Continue with workflow execution? (yes/no): " confirm
if [[ "$confirm" != "yes" ]]; then
    echo "Execution cancelled"
    exit 0
fi

echo ""
echo -e "${GREEN}=== Starting wf5 workflow ===${NC}"
echo ""

mkdir -p "$OUTDIR"

nextflow run main.nf \
    -profile azure \
    --samplesheet "$SAMPLESHEET" \
    --data_dir "$DATA_DIR" \
    --outdir "$OUTDIR" \
    --kingdom_refs "$KINGDOM_REFS" \
    --blast_db "$BLAST_DB" \
    --locus_panel "$LOCUS_PANEL" \
    ${POLISH_FLAG:-} \
    $RESUME

exit_code=$?

echo ""
if [[ $exit_code -eq 0 ]]; then
    echo -e "${GREEN}=== Workflow completed successfully ===${NC}"
    echo "Output directory: $OUTDIR"
else
    echo -e "${RED}=== Workflow failed ===${NC}"
    echo "Exit code: $exit_code"
    echo "Check .nextflow.log for details"
fi

exit $exit_code
