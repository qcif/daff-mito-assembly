# Storing data on Azure blob

Our storage accounts are used to provision data for multiple workflows in 
addition to this one. Therefore, it's important that we take care to create
blob paths that separate each workflow's resources appropriately.

Currently we are aiming for (top-level is bucket name; starred is for this workflow):

```
cache/
├── taxodactyl/
└── *wf5/

integration-fixtures/
└── *wf5/

refdata-wf4/

*refdata-wf5/

scripts/
├── taxodactyl/
└── *wf5/

workdata/
├── taxodactyl/
└── *wf5/
```
