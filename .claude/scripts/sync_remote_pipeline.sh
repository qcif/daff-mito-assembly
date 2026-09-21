#!/usr/bin/env bash
# Sync the pipeline checkout to the daff-admin remote VM.
#
# Mirrors every path `main.nf` / `nextflow.config` / `conf/*.config`
# resolve against ${projectDir} onto /mnt/data/wf5/pipeline/. Getting
# this list wrong is the classic remote-run failure: a missing
# scripts/report/ renders every per-sample report.html as a Jinja
# TemplateNotFound traceback, and a missing subworkflows/ aborts at
# include-time (both hit during task 49).
#
# Usage:
#   bash .claude/scripts/sync_remote_pipeline.sh [--with-stub-fixtures]
#
# --with-stub-fixtures also syncs tests/stub/ (only needed to run
# `-profile stub` remotely as a wiring sanity check).

set -euo pipefail

REMOTE="${WF5_REMOTE:-daff-admin}"
REMOTE_ROOT="${WF5_REMOTE_ROOT:-/mnt/data/wf5}"
REMOTE_PIPELINE="$REMOTE_ROOT/pipeline"

cd "$(git rev-parse --show-toplevel)"

WITH_STUB=0
[[ "${1:-}" == "--with-stub-fixtures" ]] && WITH_STUB=1

# Top-level files. No trailing-slash subtlety applies to plain files.
FILES=(main.nf nextflow.config)

# Directories. NOTE: no trailing slash — `rsync -az modules dest/`
# creates dest/modules/, whereas `rsync -az modules/ dest/` dumps the
# *contents* into dest/ and flattens the layout. The flattened form
# silently produces a pipeline dir that fails at include-time.
#   bin/         — all bin/*.py staged onto every process
#   conf/        — base/containers/stub/integration/azure configs
#   modules/     — process definitions
#   subworkflows/— annotation_scoring
#   assets/      — loci.json, gene sets, JSON schemas (nextflow.config)
#   scripts/     — report/templates + report/static (main.nf)
DIRS=(bin conf modules subworkflows assets scripts)

echo "Syncing pipeline -> $REMOTE:$REMOTE_PIPELINE"
ssh "$REMOTE" "mkdir -p '$REMOTE_PIPELINE'"

rsync -az "${FILES[@]}" "$REMOTE:$REMOTE_PIPELINE/"

for d in "${DIRS[@]}"; do
    # --delete keeps the remote a true mirror so an edited/removed file
    # locally can't leave a stale copy running remotely. Safe here
    # because $REMOTE_PIPELINE is a pure sync target, never a place
    # run artifacts are written — work/ and output/ live beside it
    # under $REMOTE_ROOT, and the resource-override config is synced
    # to $REMOTE_ROOT (not inside the mirror) by design.
    rsync -az --delete \
        --exclude='__pycache__/' --exclude='*.pyc' \
        "$d" "$REMOTE:$REMOTE_PIPELINE/"
    echo "  synced $d/"
done

if [[ "$WITH_STUB" == "1" ]]; then
    ssh "$REMOTE" "mkdir -p '$REMOTE_PIPELINE/tests'"
    rsync -az --delete tests/stub "$REMOTE:$REMOTE_PIPELINE/tests/"
    echo "  synced tests/stub/ (stub fixtures)"
fi

# Resource override lives OUTSIDE the mirrored pipeline dir so
# --delete above can never remove it.
if [[ -f tests/manual/remote_resources.config ]]; then
    rsync -az tests/manual/remote_resources.config \
        "$REMOTE:$REMOTE_ROOT/remote_resources.config"
    echo "  synced remote_resources.config -> $REMOTE_ROOT/"
fi

echo "Done. Pipeline synced to $REMOTE:$REMOTE_PIPELINE"
