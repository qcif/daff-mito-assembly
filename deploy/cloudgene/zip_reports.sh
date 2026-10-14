#!/usr/bin/env bash
# Cloudgene post-run step: bundle run-report.html, run_manifest.json and
# each sample's report.html into <outdir>/reports.zip, keeping the
# <sample_id>/report.html layout so the run report's links still work.

set -euo pipefail

usage() {
    echo "Usage: $0 --dir <outdir>"
}

DIR=""
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --dir)
            [[ -n "${2:-}" ]] || { echo "Missing value for --dir"; usage; exit 1; }
            DIR="${2%/}"
            shift 2
            ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
done
[[ -n "$DIR" ]] || { usage; exit 1; }

if [[ ! -f "${DIR}/run-report.html" ]]; then
    echo "ERROR: ${DIR}/run-report.html not found."
    exit 1
fi

staging="$(mktemp -d)"
trap 'rm -rf "$staging"' EXIT
mkdir "${staging}/reports"

cp "${DIR}/run-report.html" "${staging}/reports/"
[[ ! -f "${DIR}/run_manifest.json" ]] \
    || cp "${DIR}/run_manifest.json" "${staging}/reports/"

for report in "${DIR}"/*/report.html; do
    [[ -f "$report" ]] || continue
    sample_dir="$(basename "$(dirname "$report")")"
    mkdir "${staging}/reports/${sample_dir}"
    cp "$report" "${staging}/reports/${sample_dir}/"
done

rm -f "${DIR}/reports.zip"
(cd "$staging" && zip -qrX "$(realpath "$DIR")/reports.zip" reports)

echo "All workflow reports have been zipped"
echo 'Download all reports from reports.zip under the "Results" tab.'
