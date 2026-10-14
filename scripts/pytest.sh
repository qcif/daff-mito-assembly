#!/usr/bin/env bash
# Runs bin/*.py unit tests (pytest, with branch coverage) inside the
# neoformit/daff-wf5-scripts:test image, so the runtime environment
# matches production (pinned biopython etc.) rather than the host.
#
# Usage:
#   bash scripts/pytest.sh [pytest args...]
#
# Examples:
#   bash scripts/pytest.sh scripts/tests/test_plastid_canonicalise.py
#   bash scripts/pytest.sh scripts/tests/ -v
#
# Defaults to all of scripts/tests/ if no args are given. Prints a
# branch-coverage report for the modules imported by the selected
# tests.
#
# Environment:
#   COVERAGE_FAIL_UNDER  If set (e.g. 100), the coverage report fails
#                        when total coverage is below this percentage.
#                        Unset by default, so ad-hoc runs of a single
#                        test file report only. CI sets it to 100.

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
IMAGE="neoformit/daff-wf5-scripts:test"

if [[ "$(docker images -q "$IMAGE" 2>/dev/null)" == "" ]]; then
    echo "Building ${IMAGE} (scripts/Dockerfile, --build-arg TEST=1)..."
    docker build --build-arg TEST=1 -t "$IMAGE" -f "${REPO_ROOT}/scripts/Dockerfile" "${REPO_ROOT}/scripts/"
fi

TARGETS=("$@")
if [[ ${#TARGETS[@]} -eq 0 ]]; then
    TARGETS=("scripts/tests/")
fi

COVERAGE_FILE="/tmp/.coverage_$$"
FAIL_UNDER_ARG="${COVERAGE_FAIL_UNDER:+--fail-under=${COVERAGE_FAIL_UNDER}}"

docker run --rm --user "$(id -u):$(id -g)" \
    -v "${REPO_ROOT}:/repo:ro" \
    -v /tmp:/tmp \
    -e HOME=/tmp \
    -e COVERAGE_FILE="$COVERAGE_FILE" \
    -w /repo \
    "$IMAGE" \
    bash -c "
        python -m coverage run --branch -m pytest ${TARGETS[*]} -p no:cacheprovider &&
        python -m coverage report -m --include='*/bin/*.py,*/bin/**/*.py' ${FAIL_UNDER_ARG}
    "
