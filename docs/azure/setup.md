# Azure environment setup

Everything needed to recreate the wf5 Azure environment from scratch.
The environment is small: one resource group, two storage accounts,
one Batch account, and a handful of blob containers. No Terraform or
ARM templates — the resources are simple enough to provision with
`az` CLI commands documented here.

## Prerequisites

- [Azure CLI](https://learn.microsoft.com/en-us/cli/azure/install-azure-cli)
  `>= 2.60` — `az --version`
- Logged in to the correct subscription:
  ```bash
  az login
  az account set --subscription "DAFF Biosecurity"
  az account show   # confirm
  ```
- Permissions: **Contributor** on the resource group (or Owner to create it).

## 1. Resource group

```bash
az group create \
    --name daff-biosecurity \
    --location australiaeast
```

## 2. Storage accounts

Two accounts:

| Account | SKU | Purpose |
|---|---|---|
| `daffstandard` | Standard LRS | Scripts, refdata, integration fixtures, Nextflow work logs |
| `daffpremium` | Premium LRS | (Reserved for high-IOPS workloads; currently unused) |

```bash
az storage account create \
    --name daffstandard \
    --resource-group daff-biosecurity \
    --location australiaeast \
    --sku Standard_LRS \
    --kind StorageV2 \
    --https-only true \
    --min-tls-version TLS1_2

az storage account create \
    --name daffpremium \
    --resource-group daff-biosecurity \
    --location australiaeast \
    --sku Premium_LRS \
    --kind BlockBlobStorage \
    --https-only true \
    --min-tls-version TLS1_2
```

### 2.1 Blob containers (daffstandard)

| Container | Public access | Purpose |
|---|---|---|
| `scripts` | off | Batch node start-task scripts (wf5: `wf5/setup.sh`) |
| `workdata` | off | Nextflow work directories for Azure Batch runs, wf5 under `workdata/wf5/<run id>/` (14-day lifecycle) |
| `cache` | off | Nextflow cloud cache, wf5 under `cache/wf5/test/` and `cache/wf5/prod/` |
| `refdata-wf5` | blob | Versioned reference bundle: `<version>/refs.tar.gz` plus the extracted tree (`scripts/publish_refs_extracted.sh`) |
| `test-data` | blob | Tier 2 pre-recruited FASTQ fixtures for nightly CI, under `wf5/` |

```bash
# Private containers
for container in scripts workdata cache; do
    az storage container create \
        --account-name daffstandard \
        --name "$container" \
        --auth-mode login \
        --public-access off
done

# Public blob access — data is public; no credentials needed for CI downloads
for container in refdata-wf5 test-data; do
    az storage container create \
        --account-name daffstandard \
        --name "$container" \
        --auth-mode login \
        --public-access blob
done
```

The `wf5/` prefixes in `workdata` and `cache` are virtual: they appear in
listings once Nextflow writes its first blob. `run-wf5-azure.sh` checks the
`workdata` container is reachable before each run.

### 2.2 Lifecycle policy — workdata auto-delete

Nextflow work directories accumulate on blob storage. The policy in
[`deploy/azure/storage-policy.json`](../../deploy/azure/storage-policy.json)
deletes blobs under `workdata/` 14 days after last modification. It is
account-wide: other workflows share the container. Because `-resume` cache
entries point at work-dir files, this caps how far back `-resume` can reach
at 14 days. `cache/` has no rule of its own.

```bash
az storage account management-policy create \
    --account-name daffstandard \
    --resource-group daff-biosecurity \
    --policy @deploy/azure/storage-policy.json
```

## 3. Azure Batch account

```bash
az batch account create \
    --name daffwf5batch \
    --resource-group daff-biosecurity \
    --location australiaeast \
    --storage-account daffstandard
```

Note the endpoint returned (`https://daffwf5batch.australiaeast.batch.azure.com`)
— it goes into `.env.azure` as `AZURE_BATCH_ENDPOINT`.

### 3.1 Batch pools

Two pools, both created from rendered templates (fill placeholders first,
see §5):

| Pool | Template | VM |
|---|---|---|
| `wf5_test` | `pool-setup-test.json.template` | `Standard_L8as_v3` (8 vCPU / 64 GiB) |
| `wf5` | `pool-setup.json.template` | `Standard_L8as_v3` (8 vCPU / 64 GiB) |

Use `wf5_test` for development, stub and integration runs and `wf5` for
production. Create `wf5_test` first and verify against it. The VM size
must clear `process_high` (8 CPU / 32 GB) with margin; `taskSlotsPerNode`
equals the vCPU count. Image: `microsoft-dsvm:ubuntu-hpc:2404` with node
agent `batch.node.ubuntu 24.04` (check with
`az batch pool supported-images list`).

```bash
source deploy/azure/batch-helpers.sh

az_pool_create --json deploy/azure/pool-setup-test.json.ignore --autoscale
az_pool_create --json deploy/azure/pool-setup.json.ignore --autoscale
```

The start task downloads `scripts/wf5/setup.sh` and installs a pinned
`azcopy` into `$AZ_BATCH_NODE_SHARED_DIR/bin`. Reference data is not
pre-staged: Nextflow stages it per task from `az://refdata-wf5/`. Nodes
start at 0 and scale up only when tasks are queued (autoscale formula in
[`batch-helpers.sh`](../../deploy/azure/batch-helpers.sh)).

## 4. Environment file (`.env.azure`)

Create `.env.azure` at the repo root (git-ignored):

```bash
# .env.azure — fill from Azure portal / CLI output
AZURE_BATCH_ACCOUNT_NAME=daffwf5batch
AZURE_BATCH_ACCESS_KEY=<from: az batch account keys list --name daffwf5batch ...>
AZURE_BATCH_ENDPOINT=https://daffwf5batch.australiaeast.batch.azure.com
AZURE_STORAGE_ACCOUNT_KEY=<from: az storage account keys list --account-name daffstandard ...>
STORAGE_ACCOUNT_STD=daffstandard
STORAGE_CONTAINER_SCRIPTS=scripts
ANALYST_NAME="<your name>"
FACILITY_NAME="QCIF"
```

Retrieve the keys:

```bash
az batch account keys list \
    --name daffwf5batch \
    --resource-group daff-biosecurity \
    --query "primary" -o tsv

az storage account keys list \
    --account-name daffstandard \
    --query "[0].value" -o tsv
```

## 5. Start-task script and SAS token

The `test-data` and `refdata-wf5` containers are public, so neither CI nor
Batch tasks need a token for them. The only token is a short-lived read
token so a pool can download its start-task script from the private
`scripts` container.

Upload the script (it has no placeholders), then render both pool
templates with a 30-day token:

```bash
source deploy/azure/batch-helpers.sh
az_storage_upload deploy/azure/setup.sh wf5/setup.sh

az_sas_generate wf5/setup.sh daffstandard scripts 30
# paste the token over <SETUP_SCRIPT_SAS_TOKEN> in each template, saving as
#   deploy/azure/pool-setup.json.ignore
#   deploy/azure/pool-setup-test.json.ignore
```

## 6. GitHub secrets

No GitHub secrets are required for CI. Both `test-data` and
`refdata-wf5` use public blob access; `integration.yml` fetches them with plain
`curl` — no SAS tokens.

Confirm: `gh secret list`

## 7. Token renewal schedule

| Token | Expiry | Where stored |
|---|---|---|
| Pool JSON setup-script token | 30 days from creation | `pool-setup*.json.ignore` (used at pool creation and node start) |

Nodes started after the token expires cannot download the start-task
script. To renew: re-run §5 and update the start task in place
(`az_pool_update --json`) or recreate the pool.

## 8. Verify the environment

```bash
# Storage containers exist
az storage container list \
    --account-name daffstandard \
    --auth-mode login \
    --query "[].name" -o tsv

# Batch pool state
source deploy/azure/batch-helpers.sh
az_pool_show

# Integration fixtures present
az storage blob list \
    --account-name daffstandard \
    --container-name test-data \
    --prefix wf5/v2026.07/ \
    --auth-mode login \
    --query "[].{name:name, size:properties.contentLength}" \
    -o table
```
