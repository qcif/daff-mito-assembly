#!/usr/bin/env bash
#
# Azure Batch start task for the wf5 organelle assembly pipeline.
#
# Installs a pinned azcopy into $AZ_BATCH_NODE_SHARED_DIR/bin. Nextflow's
# Azure Batch task wrapper calls azcopy to stage each task's inputs and
# outputs, and conf/azure.config sets copyToolInstallMode = 'off', so the
# pool must provide it. That directory is mounted into task containers.
#
# The real binary is installed as azcopy.real behind a wrapper named azcopy.
# Task users have no writable home directory on this image, and azcopy
# otherwise fails creating ~/.azcopy, so the wrapper points its job plan
# and log dirs at /tmp.
#
# azcopy version: 10.27.1
#
# No placeholders or secrets: upload as-is to scripts/wf5/setup.sh.

set -x
set -euo pipefail

AZCOPY_VERSION="10.27.1"
AZCOPY_URL="https://github.com/Azure/azure-storage-azcopy/releases/download/v${AZCOPY_VERSION}/azcopy_linux_amd64_${AZCOPY_VERSION}.tar.gz"
INSTALL_DIR="${AZ_BATCH_NODE_SHARED_DIR:?AZ_BATCH_NODE_SHARED_DIR not set}/bin"

install_azcopy() {
    if [[ -x "${INSTALL_DIR}/azcopy.real" ]] \
        && "${INSTALL_DIR}/azcopy.real" --version | grep -q "${AZCOPY_VERSION}"; then
        echo "azcopy ${AZCOPY_VERSION} already installed" >&2
    else
        local tmp
        tmp="$(mktemp -d)"
        wget -q -O "${tmp}/azcopy.tar.gz" "${AZCOPY_URL}"
        tar -xf "${tmp}/azcopy.tar.gz" -C "${tmp}" --strip-components=1
        mkdir -p "${INSTALL_DIR}"
        install -m 0755 "${tmp}/azcopy" "${INSTALL_DIR}/azcopy.real"
        rm -rf "${tmp}"
    fi

    cat > "${INSTALL_DIR}/azcopy" <<'WRAPPER'
#!/usr/bin/env bash
export AZCOPY_JOB_PLAN_LOCATION="/tmp/azcopy-$(id -u)/plans"
export AZCOPY_LOG_LOCATION="/tmp/azcopy-$(id -u)/logs"
mkdir -p "${AZCOPY_JOB_PLAN_LOCATION}" "${AZCOPY_LOG_LOCATION}"
exec "$(dirname "$0")/azcopy.real" "$@"
WRAPPER
    chmod 0755 "${INSTALL_DIR}/azcopy"
    "${INSTALL_DIR}/azcopy" --version >&2
}

echo "===== START TASK: installing azcopy ${AZCOPY_VERSION} =====" >&2
install_azcopy
echo "===== START TASK: completed =====" >&2
