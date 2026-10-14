# Storing data on Azure blob

Our storage accounts are used to provision data for multiple workflows in 
addition to this one. Therefore, it's important that we take care to create
blob paths that separate each workflow's resources appropriately.

The below are using the `daffstandard` storage account (standard blob).

Currently we are aiming for (top-level is bucket name; starred is for this workflow):

```
cache/
├── taxodactyl/
└── *wf5/

refdata-wf4/

*refdata-wf5/

scripts/
├── taxodactyl/
└── *wf5/

test-data/
├── taxodactyl/
└── *wf5/

workdata/
├── taxodactyl/
└── *wf5/
```

## Prefixes used by this workflow

| Prefix | Used by |
|---|---|
| `refdata-wf5/<version>/` | Tarball plus extracted bundle (spec 04 §4.4–4.6); staged per task on Azure Batch |
| `scripts/wf5/setup.sh` | Start task for the `wf5` and `wf5_test` Batch pools (installs `azcopy`) |
| `workdata/wf5/<run id>/` | Nextflow work dir, one per run (`deploy/azure/run-wf5-azure.sh`); 14-day lifecycle rule |
| `cache/wf5/{test,prod}/` | Nextflow cloud cache for `-resume`, split by pool |
| `test-data/wf5/<fixtures version>/` | Integration fixtures |

Empty directories do not survive blob storage; see the "Directory outputs are never empty" rule in spec/01-pipeline-flow.md §1a.
