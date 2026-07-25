# Azure Batch deployment

The pipeline supports Azure Batch as an execution backend (`-profile azure`).
Deployment templates and helpers live under [`deploy/azure/`](../../deploy/azure/).

## Repo layout

```
deploy/azure/
├── setup.sh.template            # Start task: install azcopy + stage refdata bundle
├── pool-setup.json.template     # Production pool (D8s_v3, autoscaling to 1 node)
├── pool-setup-test.json.template  # Test pool (D4s_v3)
├── batch-helpers.sh             # Shell functions for pool/job/storage management
├── run-wf5.sh                   # Invoke the pipeline against the Azure pool
├── storage-policy.json          # Lifecycle rule: expire workdata/ blobs after 14 days
└── .gitignore                   # Rendered *.ignore files are git-ignored
```

`.template` files carry placeholder tokens (`<REFDATA_SAS_TOKEN>`,
`<SETUP_SCRIPT_SAS_TOKEN>`, `<REFDATA_BUNDLE_VERSION>`) that must be filled
before use. Rendered copies are named `*.ignore` and git-ignored.

## Prerequisites

Set the following environment variables (e.g. in `.env.azure`, git-ignored):

```bash
AZURE_BATCH_ACCOUNT_NAME=...
AZURE_BATCH_ACCESS_KEY=...
AZURE_BATCH_ENDPOINT=...
AZURE_STORAGE_ACCOUNT_KEY=...
STORAGE_ACCOUNT_STD=daffstandard
STORAGE_ACCOUNT_PREM=daffpremium
STORAGE_CONTAINER_SCRIPTS=scripts
ANALYST_NAME="..."
FACILITY_NAME="..."
```

## First-time setup

1. **Build the reference bundle** per [plan.md §4.4](../../plan.md) — outputs
   `refs/v<YYYY.MM>/`.
2. **Upload** the bundle to the `refdata-wf5` blob container on the
   `daffstandard` storage account.
3. **Render** `setup.sh.template` → `setup.sh.ignore`: fill
   `<REFDATA_SAS_TOKEN>` (generate with `az_sas_generate "" daffstandard refdata-wf5`)
   and `<REFDATA_BUNDLE_VERSION>` (e.g. `v2026.06`).
4. **Upload** `setup.sh.ignore` to `daffstandard/scripts/setup-wf5.sh`
   (`az_storage_upload deploy/azure/setup.sh.ignore setup-wf5.sh`).
5. **Render** `pool-setup.json.template` → `pool-setup.json.ignore`: fill
   `<SETUP_SCRIPT_SAS_TOKEN>` for the setup script blob
   (`az_sas_generate setup-wf5.sh`).
6. **Create the pool** with autoscaling enabled:
   `az_pool_create --json deploy/azure/pool-setup.json.ignore --autoscale`

## Running the pipeline

```bash
source deploy/azure/batch-helpers.sh   # loads .env.azure and helpers
./deploy/azure/run-wf5.sh \
    --samplesheet samples.csv \
    --data-dir /path/to/reads \
    --outdir output/run1
```

The pool auto-scales up when tasks are queued and back down after 15 minutes
of idle time.

## Debugging

- `az_pool_show` — pool state and node counts
- `az_node_logs` — start task stdout/stderr from a live node
- `az_job_logs` — task status for the most recent job
- `az_fetch_nf_logs .nextflow.log ./logs` — pull `.command.{out,err,sh}` and
  exit codes for every task in a run

## Notes

- All process I/O must flow through channels ([plan.md §1a](../../plan.md)) —
  Azure Batch only stages files it can see on a channel.
- Container images are pulled from Docker Hub
  ([`neoformit/daff-wf5-scripts`](https://hub.docker.com/r/neoformit/daff-wf5-scripts)).
- Reference bundle version is baked into the pool's start task; bumping the
  bundle version requires re-rendering `setup.sh.ignore`, re-uploading, and
  recreating the pool (or `az_pool_update --json` + `az_pool_resize 0 && 1`
  to force node recreation).
