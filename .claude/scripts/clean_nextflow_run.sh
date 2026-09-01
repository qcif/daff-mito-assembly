#!/usr/bin/env bash
# Clean up local Nextflow run artifacts (logs, work dir, stub output)
# left behind after ad-hoc `nextflow run` invocations during task work.
#
# Per CLAUDE.md: keep only the two most recent .nextflow.log.N files,
# and remove ./work/ entirely when a clean run is wanted rather than a
# `-resume`.
#
# Usage: bash .claude/scripts/clean_nextflow_run.sh

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

# Keep .nextflow.log and .nextflow.log.1/.2, drop anything older.
for f in .nextflow.log.[3-9] .nextflow.log.[1-9][0-9]*; do
    [[ -e "$f" ]] && rm -f "$f"
done

rm -rf work
rm -rf tests/output

echo "Cleaned: work/, tests/output/, stale .nextflow.log.N files."
