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
| `scripts` | off | Batch node start-task scripts (e.g. `setup-wf5.sh`) |
| `workdata` | off | Nextflow work directory for Azure Batch runs (14-day lifecycle) |
| `cache` | off | Nextflow cache for Azure Batch runs |
| `refdata-wf5` | blob | Versioned production reference bundle (published by `scripts/build_refs.sh`) |
| `integration-fixtures` | blob | Tier 2 pre-recruited FASTQ fixtures for nightly CI |

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
for container in refdata-wf5 integration-fixtures; do
    az storage container create \
        --account-name daffstandard \
        --name "$container" \
        --auth-mode login \
        --public-access blob
done
```

### 2.2 Lifecycle policy — workdata auto-delete

Nextflow work directories accumulate on blob storage. The policy in
[`deploy/azure/storage-policy.json`](../../deploy/azure/storage-policy.json)
deletes blobs under `workdata/` 14 days after last modification.

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

### 3.1 Batch pool

The pool is created from a rendered template. Fill placeholders
first (see §5), then:

```bash
source deploy/azure/batch-helpers.sh

# Production pool (D8s_v3, auto-scales to 1 node on demand)
az_pool_create --json deploy/azure/pool-setup.json.ignore --autoscale

# Optional: test pool (D4s_v3)
az_pool_create --json deploy/azure/pool-setup-test.json.ignore
```

The start task downloads `setup-wf5.sh` from blob storage and stages
the refdata bundle to `/mnt/refdata` on each node. Nodes start at 0
and scale up only when tasks are queued (autoscale formula in
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
STORAGE_ACCOUNT_PREM=daffpremium
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

## 5. SAS tokens

SAS tokens are only needed for Batch node operations. The `integration-fixtures`
and `refdata-wf5` containers are public (`--public-access blob`), so CI fetches
require no credentials.

### 5.1 Refdata token (Batch start task)

The pool start task (`setup-wf5.sh`) needs a token to download the
refdata bundle onto each node. This is a blob-level read token with
a longer expiry (the pool's lifetime):

```bash
source deploy/azure/batch-helpers.sh
az_sas_generate "" daffstandard refdata-wf5 730   # 2-year expiry
```

Paste the token into `deploy/azure/setup.sh.template` as
`<REFDATA_SAS_TOKEN>`, save as `setup.sh.ignore`, upload:

```bash
az storage blob upload \
    --account-name daffstandard \
    --container-name scripts \
    --name setup-wf5.sh \
    --file deploy/azure/setup.sh.ignore \
    --auth-mode login \
    --overwrite
```

### 5.4 Setup-script token (pool JSON)

The pool JSON needs a short-lived read token to download the start
task script during pool creation:

```bash
az_sas_generate setup-wf5.sh daffstandard scripts 30   # 30-day expiry
```

Paste into `pool-setup.json.template` as `<SETUP_SCRIPT_SAS_TOKEN>`,
save as `pool-setup.json.ignore`.

## 6. GitHub secrets

No GitHub secrets are required for CI. Both `integration-fixtures` and
`refdata-wf5` use public blob access; `integration.yml` fetches them with plain
`curl` — no SAS tokens.

Confirm: `gh secret list`

## 7. Token renewal schedule

SAS tokens expire. Calendar reminders recommended:

| Token | Expiry | Where stored |
|---|---|---|
| Batch start-task refdata token | 2 years from creation | `setup.sh.ignore` → blob |
| Pool JSON setup-script token | 30 days from creation | `pool-setup.json.ignore` (one-off) |
Renew the start-task token: re-run §5.3, re-upload `setup.sh.ignore`,
and recreate the pool (or update the start task in place with
`az batch pool set`).

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
    --container-name integration-fixtures \
    --prefix wf5/v2026.07/ \
    --auth-mode login \
    --query "[].{name:name, size:properties.contentLength}" \
    -o table
```
