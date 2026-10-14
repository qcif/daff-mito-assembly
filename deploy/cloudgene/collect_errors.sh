#!/usr/bin/env bash
# Cloudgene post-run step: summarise sample status, print the tail of
# each failed task's logs and zip <outdir>/errors/ to errors.zip.

set -euo pipefail

EXIT_NO_MANIFEST=3
TAIL_LINES=10

usage() {
    echo "Usage: $0 --dir <outdir>"
    echo "Reads <outdir>/run_manifest.json and <outdir>/errors/"
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

manifest="${DIR}/run_manifest.json"
errors_dir="${DIR}/errors"
errors_zip="${DIR}/errors.zip"

if [[ ! -s "$manifest" ]]; then
    echo "ERROR: ${manifest} not found or empty."
    echo "The run did not finish writing its results; see the workflow log."
    exit "$EXIT_NO_MANIFEST"
fi

# Status counts, then each error/invalid_input sample with its reason.
python3 - "$manifest" <<'PY'
import json
import sys

manifest = json.load(open(sys.argv[1]))
print('Sample status summary')
for status, count in manifest['status_breakdown'].items():
    if count:
        print(f'  {status}: {count}')
for sample in manifest['samples']:
    if sample['sample_status'] in ('error', 'invalid_input'):
        print(f"\n{sample['sample_id']} [{sample['sample_status']}]: "
              f"{sample.get('sample_status_reason')}")
PY

if [[ ! -d "$errors_dir" ]]; then
    echo
    echo "No task errors"
    rm -f "$errors_zip"
    exit 0
fi

# Logs grouped by sample: <PROCESS>/<sample_id>.<ext>
echo
echo "Task logs (last ${TAIL_LINES} lines)"
while IFS= read -r sample; do
    echo
    echo "=== ${sample} ==="
    while IFS= read -r -d '' log_file; do
        echo "--- ${log_file#"$errors_dir"/}"
        tail -n "$TAIL_LINES" "$log_file" || true
    done < <(find "$errors_dir" -type f \
        \( -name "${sample}.err" -o -name "${sample}.log" \) \
        -print0 | sort -z)
done < <(find "$errors_dir" -type f \( -name '*.err' -o -name '*.log' \) \
    -printf '%f\n' | sed -E 's/\.(err|log)$//' | sort -u)

rm -f "$errors_zip"
(cd "$DIR" && zip -qr errors.zip errors)

echo
echo "Zipped task logs to: ${errors_zip}"
echo 'Download all task logs from errors.zip under the "Results" tab.'
