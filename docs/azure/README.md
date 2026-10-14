# Azure Batch deployment

The pipeline runs on Azure Batch with `-profile azure`. Nextflow runs on a
small launcher machine; each process becomes a Batch task on a named pool.
Every input and output moves through Azure Blob Storage.

**New to this environment?** See [setup.md](setup.md) for from-scratch
provisioning (resource group, storage, Batch account, containers, pools,
SAS tokens).

## Repo layout

```
deploy/azure/
├── setup.sh                     # Pool start task: install pinned azcopy only
├── pool-setup.json.template     # Production pool `wf5` (L8as_v3)
├── pool-setup-test.json.template  # Test pool `wf5_test` (L8as_v3)
├── batch-helpers.sh             # Shell functions for pool/job/storage management
├── run-wf5-azure.sh             # Launcher: stub, integration and sample runs
└── storage-policy.json          # Lifecycle rule: expire workdata/ blobs after 14 days
conf/azure.config                # Nextflow profile (named pool, az:// work dir + cache)
scripts/publish_refs_extracted.sh  # Publish an unpacked refdata bundle to blob
```

`.template` files carry a `<SETUP_SCRIPT_SAS_TOKEN>` placeholder. Rendered
copies are named `*.ignore` and are git-ignored.

## Blob layout

All on `daffstandard` (see spec/04a-azure-blob-storage.md):

| Location | Used for |
|---|---|
| `refdata-wf5/<version>/refs.tar.gz` | Tarball for `scripts/fetch_refs.sh` and CI |
| `refdata-wf5/<version>/{recruit,validate,proteins,annotate,manifest.json}` | Extracted bundle; Nextflow stages it into each task |
| `scripts/wf5/setup.sh` | Pool start task |
| `workdata/wf5/<run id>/` | Nextflow work dir, one per run. Expires 14 days after last write, so `-resume` cannot reach further back than that |
| `cache/wf5/{test,prod}/` | Nextflow cloud cache (`-resume` state). Entries are only useful while their work dir exists. No lifecycle rule: entries are tiny, and stale ones are harmless |
| `test-data/wf5/<version>/` | Integration fixtures |

Reference data is staged **per task** from `az://refdata-wf5/<version>/`,
not pre-staged on nodes. The largest directories are `recruit/` (~955 MB)
and `validate/` (~936 MB).

## Pools

| Pool | VM | Use |
|---|---|---|
| `wf5_test` | `Standard_L8as_v3`, 8 vCPU / 64 GiB | Development, stub and integration runs. Default |
| `wf5` | `Standard_L8as_v3`, 8 vCPU / 64 GiB | Production runs. Must be asked for with `--pool wf5` |

Both fit a `process_high` task (8 CPU / 32 GB) with headroom, and
`taskSlotsPerNode` equals the vCPU count so Batch packing matches
Nextflow's `cpus`. They differ only in id. Both draw on the account's 80-core LASv3 quota,
shared with other workflows' pools. Image:
`microsoft-dsvm:ubuntu-hpc:2404` (Ubuntu 24.04, Docker preinstalled).

The start task installs a pinned `azcopy` into
`$AZ_BATCH_NODE_SHARED_DIR/bin`. Nextflow's Azure executor calls it to stage
task inputs and outputs. `conf/azure.config` sets
`copyToolInstallMode = 'off'`, so the pool is responsible for providing it.

## Prerequisites

`.env.azure` at the repo root (git-ignored), see `.env.sample`. The
variables used by Nextflow and the launcher:

```bash
AZURE_BATCH_ACCOUNT_NAME=...
AZURE_BATCH_ENDPOINT=https://<account>.<region>.batch.azure.com
AZURE_BATCH_ACCESS_KEY=...
STORAGE_ACCOUNT_STD=daffstandard
AZURE_STORAGE_ACCOUNT_KEY=...
```

## Running the pipeline

```bash
# Cheap wiring / auth / pool check
deploy/azure/run-wf5-azure.sh --stub --refs-version v2026.09_1

# Integration fixtures; runs tests/integration/assertions.sh on success
deploy/azure/run-wf5-azure.sh --integration --refs-version v2026.09_1

# Your own samples, local reads (uploaded once into the work dir)
deploy/azure/run-wf5-azure.sh --refs-version v2026.09_1 \
    --samplesheet samples.csv --data-dir /path/to/reads --outdir out

# Samplesheet listing az:// read URIs: no --data-dir, nothing uploaded
deploy/azure/run-wf5-azure.sh --refs-version v2026.09_1 \
    --samplesheet samples_az.csv --outdir out
```

`--refs-version` is required. `--pool` defaults to `wf5_test`. Other flags:
`--polish`, `-resume` (see below), `--trace` (needs `ps` in the task containers, which most wf5 images lack), `--yes` (skip confirmation), `--fixtures-version`.
After Nextflow exits 0 the launcher reads `run_manifest.json` and prints a
per-status sample count. It then lists every `error` or `invalid_input`
sample with its reason and exits **3**. Exit **4** means no
`run_manifest.json` was produced. Any other non-zero code is Nextflow's own
(or a failed `--integration` assertion, which is run first). The
integration mode drops `INT-CORRUPT-01`, a locally generated fixture, from
its samplesheet.
The launcher checks the pool is `active` and the `workdata` container is
reachable before starting Nextflow, and prints the work dir and the
`az_fetch_nf_logs` command for failures.

## Resuming a run

`-resume` state is in `cache/wf5/{test,prod}/<session id>/`, so it survives
a fresh checkout, but Nextflow 25.04 still needs the session id and a local
history entry. The session id is in the original run's `.nextflow.log`
(`Session UUID`). On a new checkout:

```bash
mkdir -p .nextflow
printf '2026-01-01 00:00:00\t-\tresume\tOK\t-\t<session id>\tnextflow run main.nf\n' \
    > .nextflow/history
deploy/azure/run-wf5-azure.sh <same args> -- -resume <session id>
```

Bare `-resume` finds nothing without the local history. Tasks that ignore
errors are re-run. Cache entries are only useful while the original work dir
still exists (14-day lifecycle).

## Publishing a reference bundle

After building and publishing the tarball, publish the extracted tree:

```bash
bash scripts/fetch_refs.sh v2026.09_1          # if refs/v2026.09_1 is absent
bash scripts/publish_refs_extracted.sh v2026.09_1
```

The script refuses to overwrite an existing version, uploads `manifest.json`
last, then checks the blob listing against the manifest's artefacts.

## Debugging

```bash
source deploy/azure/batch-helpers.sh
az_pool_show                                  # pool state and node counts
az_node_logs                                  # start task stdout/stderr
az_job_logs                                   # task status for the latest job
az_fetch_nf_logs .nextflow.log ./logs         # .command.* and exit codes per task
```

## Notes

- All process I/O must flow through channels (spec/01-pipeline-flow.md §1a);
  Batch only stages files it can see on a channel.
- Container images are pulled on demand from Docker Hub and Quay.
