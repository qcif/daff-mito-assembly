# Task 57 — Run the workflow on Azure Batch

**Phase:** P5 — production deployment (spec/06-phases.md).
**Goal:** Make `-profile azure` work end to end against a named Batch
pool created from `deploy/azure/pool-setup.json.template`. Add a
launcher script, `deploy/azure/run-wf5-azure.sh`, that runs the stub
and integration fixtures through that pool, plus any user samplesheet.

**Depends on:** task 56_samplesheet_read_sources.md. Without it,
read paths resolved inside `VALIDATE_SAMPLESHEET` are node-local and
the run cannot stage reads. Task 56 also defines the two read-source
modes this task's launcher must support: a local `--data_dir` (mode A)
and `az://` URIs in the sheet (mode B).

**Related spec / docs:** spec/04a-azure-blob-storage.md (blob layout),
spec/04-reference-data.md §4.4–4.5 (bundle build and fetch),
spec/01-pipeline-flow.md §1a (staging), spec/05-test-data.md (fixtures),
docs/azure/README.md, docs/azure/setup.md. CONSTITUTION rules 10, 11,
13, 15, 17, 18.

---

## 1. Overview

So far the pipeline has run on developer machines, on the `daff-admin`
VM and in GitHub Actions, all with Docker on one host. Production is
meant to run on **Azure Batch**. There, Nextflow runs on a small
launcher machine. Each process becomes a Batch task on a pool of VMs.
Every input and output moves through Azure Blob Storage, because the
VMs share no filesystem.

The repo already has Azure scaffolding: `deploy/azure/`,
`conf/azure.config` and `docs/azure/`. It was copied from an earlier
workflow and was never run against wf5. It predates most of the
pipeline. A read-through finds it would not work as written:

| Piece | Problem |
|---|---|
| `conf/azure.config` | Uses an auto-pool (`autoPoolMode = true`, D4s_v3 × 5–20) instead of the templated named pool. Sets no `workDir`. Reads env vars (`AZURE_BATCH_ACCOUNT_KEY`, `AZURE_STORAGE_ACCOUNT_NAME`, `AZURE_STORAGE_ACCESS_KEY`) that `.env.azure`, `docs/azure` and `run-wf5.sh` do not define. Those files use `AZURE_BATCH_ACCESS_KEY`, `STORAGE_ACCOUNT_STD` and `AZURE_STORAGE_ACCOUNT_KEY`. |
| `deploy/azure/run-wf5.sh` | Passes `--kingdom_refs` and `--locus_panel` (retired params) and never sets `organelle_refs`, `protein_panel`, `annotate_refs` or `refs_manifest`, so `validateParams()` aborts. Hard-codes `v2026.06`. References `$POLISH_FLAG` under `set -u` with no default. |
| `pool-setup*.json.template` | Ubuntu 20.04 `ubuntu-server-container` image and `batch.node.ubuntu 20.04` agent, which Azure Batch has retired. |
| `setup.sh.template` | Start task that installs `azcopy` to `/usr/local/bin` and pre-stages the refdata bundle to `/mnt/refdata` on each node. Nextflow cannot use a node-local path as a stageable input (§3.2), so the staged copy would go unused. Its SAS token is also unnecessary, because `refdata-wf5` is a public container. `/usr/local/bin` is likely not visible inside task containers (§3.2). |
| Blob layout | spec 04a reserves `workdata/wf5/` and `cache/wf5/`. None of the scaffolding uses those prefixes. |
| `refdata-wf5` | Holds only `<version>/refs.tar.gz`. Nextflow cannot stage a directory from inside a tarball. |

This task fixes those pieces and proves the result: a stub run and then
the three integration fixtures run on Azure Batch, and the integration
assertions pass against the output.

## 2. Scope

**In scope:** pool templates (prod + test), cutting the start task
down to an `azcopy` install, `conf/azure.config`,
publishing an extracted refdata bundle to blob, the
`run-wf5-azure.sh` launcher, the blob prefixes the run needs, and doc
and spec updates.

**Out of scope:**

- Running Azure in CI. Rule 15 allows three test surfaces. This task
  adds an on-demand manual run, not a fourth surface.
- Autoscale-formula tuning, cost reporting and spot/low-priority nodes.
- Production sample runs. A follow-up task will use `run-wf5-azure.sh`
  for those.
- Taxodactyl hand-off from blob.

## 3. Design decisions

### 3.1 Named pool from the template, not an auto-pool

Create the pool once with `az_pool_create` from a rendered
`pool-setup.json.ignore`. Point Nextflow at it by pool id:
`process.queue = '<pool id>'`, `autoPoolMode = false` and
`allowPoolCreation = false`. A named pool keeps the VM image and
sizing under version control, in one reviewed template. It
also makes "which machine ran this" answerable from the template's git
history (rule 18). With an auto-pool, those settings live inside
Nextflow config defaults.

Template changes:

- **Image:** move to a currently supported Ubuntu 22.04 (or newer)
  image with Docker preinstalled, plus the matching
  `batch.node.ubuntu 22.04` agent. Check the choice against
  `az batch pool supported-images list` at implementation time, and
  record the exact publisher/offer/sku in the template. Keep
  `containerConfiguration.type = dockerCompatible`. Pre-pull nothing
  beyond what the image needs: the biocontainer set is large and is
  pulled on demand.
- **VM size:** `METAFLYE` (and `MEDAKA` when `--polish`) request
  `process_high` = 8 CPU / 32 GB. A `Standard_D8s_v3` has 8 vCPU /
  32 GiB, and the OS and Batch agent take part of that, so a full-size
  `process_high` task has no headroom. Pick a size whose usable memory
  clears `process_high` with margin, such as an E-series 8 vCPU / 64 GiB.
  Alternatively, add `process.resourceLimits` in `conf/azure.config`
  matching the chosen VM, so Nextflow clamps requests rather than
  scheduling a task the node cannot hold. Do whichever matches the chosen
  VM, and note the reasoning in the template header comment.
  `taskSlotsPerNode` must equal the VM's vCPU count so slot-based
  packing matches Nextflow's `cpus`.
- **`plant_mt` memory:** fold in the tasks/todo.md item from task 28.
  `plant_mt.mmi` is 536 MB, and `RECRUIT`, the C2 coverage split and C3
  binning all load it. Read peak RSS from the trace of the
  `INT-PLANT-01-mt` run (§5.3) and confirm those processes fit their
  labels on the chosen VM. If they don't, raise the label for those
  processes. Record the figures in this task's outcome notes.
- **Start task installs `azcopy` only.** Keep the `startTask` block,
  pointing at the reduced script (§3.2). Remove the `/mnt/refdata`
  `containerVolumes` mount.
- **Two pools, kept separate:** `wf5` (production, from
  `pool-setup.json.template`) and `wf5_test` (development, from
  `pool-setup-test.json.template`). Both templates get the same image,
  agent and start-task changes. They differ only in id, VM size and
  `taskSlotsPerNode`. The test pool must still fit `process_high` (or
  `conf/azure.config`'s `resourceLimits` must clamp to it), or
  integration runs on it will not schedule. Dev, stub and integration
  runs go to `wf5_test`; production runs go to `wf5`. The launcher
  enforces this (§4 item 4). Both pools use the same start-task
  script.

### 3.2 Reference data is staged per task; the start task only installs `azcopy`

Decided with the user (2026-10-12): stop pre-staging reference data on
nodes and let Nextflow stage it into each task that needs it. The
bundle is small (~2 GB, with the largest single directory under 1 GB),
so a per-task copy within one region is cheap. That beats keeping
node-level setup that the pipeline can't consume anyway. The pool's
start task stays, cut down to installing `azcopy`. That tool is
infrastructure, so it belongs with the pool, not the workflow.

Why the node-local copy can't be consumed: the existing start task
copies the bundle to `/mnt/refdata` on each node, and the old run
script passed `/mnt/refdata/...` strings as params. That cannot work
under the current pipeline:

- `validateParams()` calls `file(value).exists()` on the **launcher**,
  where `/mnt/refdata` does not exist, so the run aborts.
- If that check is bypassed, `file('/mnt/refdata/...')` is a local path
  on the launcher. Nextflow tries to upload it to the work dir as a
  foreign file and fails.
- Getting around that means passing bare strings into `script:`
  blocks. Rule 13 forbids that, and the stageable-file lint fails it.

So the reference params point at **blob paths**, such as
`az://refdata-wf5/v2026.09_1/recruit`. Nextflow then stages them into
each task like any other input. That is fully rule-13 compliant and
needs no pipeline code change. Costs and consequences:

- **Per-task download.** Each task that takes a reference directory
  copies it from blob. The worst cases are `recruit/` (~955 MB) and
  `validate/` (~936 MB). That is seconds per task within one region,
  with no egress charge. It is acceptable for the expected batch sizes.
  Measure it in §5.3 and record the result. If it proves material,
  open a follow-up task. Do not reintroduce node-local strings.
- **An extracted bundle must exist on blob.** Today `refdata-wf5` holds
  only `<version>/refs.tar.gz`. Publish the unpacked tree alongside it,
  for example `refdata-wf5/<version>/{recruit,validate,proteins,annotate,manifest.json}`.
  Keep the tarball for `scripts/fetch_refs.sh` and CI. Add the
  extracted upload to whatever step publishes the tarball
  (`scripts/build_refs.sh` or its publish helper). Make it idempotent.
  Run the same completeness check `fetch_refs.sh` performs, derived
  from `manifest.json`'s `artefacts`, against the blob listing after
  upload. A partially uploaded bundle must fail loudly, as task 41's
  incident showed. Do not overwrite an existing version prefix: bundles
  are immutable (rule 10).
- **Provenance is unchanged.** `refs_manifest` becomes
  `az://refdata-wf5/<version>/manifest.json` and flows into COLLATE /
  RUN_REPORT exactly as now.

**Why `azcopy` is needed at all.** No wf5 process calls it. Nextflow's
Azure Batch executor does: each task's wrapper script uses `azcopy` to
download the task's inputs from blob into its working directory, and to
upload outputs and logs back to the `az://` work dir.

**The start task owns the install; Nextflow is told not to.** Set
`azure.batch.copyToolInstallMode = 'off'` in `conf/azure.config`. The
modes:

- `node`: Nextflow installs `azcopy` once per node, but only through
  the start task of a pool *it* creates.
- `task`: Nextflow downloads it at the start of every task.
- `off`: Nextflow assumes it is already present.

`off` makes the split of responsibility explicit. Check these
semantics against the Nextflow 25.04 docs and record what you find.

**`azcopy` must be visible inside task containers, not just on the
host.** Nextflow runs the task wrapper inside the process container.
The old script installed to `/usr/local/bin` on the host, which a
container does not see. Find out where Nextflow expects the tool.
Expect `$AZ_BATCH_NODE_SHARED_DIR/bin`: that directory is mounted into
task containers and is on the wrapper's `PATH`, and it is where
Nextflow's own `node` mode installs. Then install there. **Verify this
first**, from the generated `.command.run` / Batch task container
settings of the §5.2 stub run, before relying on it. If `off` plus a
start-task install cannot be made visible to tasks, stop and raise it
with the user rather than switching to `task` mode quietly.

Start-task changes:

- `deploy/azure/setup.sh.template` keeps only the `azcopy` install. It
  stays idempotent and installs to the location established above.
  Pin the `azcopy` version instead of the floating
  `aka.ms/downloadazcopy-v10-linux` link (rule 11), and record the
  version in the script header.
- Delete the refdata staging function, `REFDATA_*` variables, and the
  `<REFDATA_SAS_TOKEN>` / `<REFDATA_BUNDLE_VERSION>` placeholders.
  Delete the `/mnt/refdata` volume from both pool templates (§3.1).
- With no secrets left, the script needs no rendering step. Rename it
  to a plain `deploy/azure/setup.sh` if that reads more honestly.
  Upload it to `scripts/wf5/setup.sh` (spec 04a), replacing the
  root-level `scripts/setup-wf5.sh` and `scripts/setup-wf5-test.sh`.
  Check what is there first, and delete only those two. Both pools
  point at the one script.
- The pool templates still need a read SAS token for that script while
  `scripts` is a private container (`<SETUP_SCRIPT_SAS_TOKEN>`, 30-day
  expiry, used only at pool creation and node start). Keep that flow.
  Drop the 2-year refdata token row from docs/azure/setup.md §7.

The bundle version is no longer baked into the pool. It is a per-run
launcher argument (§4 item 4), so bumping it no longer means
recreating the pool.

### 3.3 Blob layout: work dir, cache, scripts

Follow spec/04a-azure-blob-storage.md exactly:

| Purpose | Location |
|---|---|
| Nextflow `workDir` | `az://workdata/wf5/` (set in `conf/azure.config`; the launcher may append a per-run suffix — see §4.4) |
| Reference bundle | `az://refdata-wf5/<version>/…` (§3.2) |
| Start-task script | `scripts/wf5/setup.sh` |
| Integration fixtures | `az://test-data/wf5/<fixtures version>/` (already present) |

The `workdata` container already exists (docs/azure/setup.md §2.1).
Blob "directories" are virtual prefixes, so `workdata/wf5/` appears
when Nextflow writes its first file. Still, make it **explicit and
checked**:

- Add the `wf5/` prefix to the provisioning steps in docs/azure/setup.md
  §2.1. If a placeholder blob is needed so the prefix shows up in
  listings and portal views, say so.
- Have `run-wf5-azure.sh` confirm in preflight that the `workdata`
  container is reachable with the configured key. Write access is
  checked implicitly on the first task, but an early, clear error is
  better.
- The existing lifecycle rule (`storage-policy.json`, prefix
  `workdata/`) already covers `workdata/wf5/`. Leave the rule as it is.
  It is account-wide and other workflows share that container.
  Confirm 14 days does not break `-resume` for a realistic re-run
  window, and note that limit in docs/azure.

**Cache: `az://cache/wf5/`.** Use Nextflow's cloud cache (the
`nf-cloudcache` plugin / `NXF_CLOUDCACHE_PATH`; check the 25.04
config syntax), so `-resume` state lives next to the work dir rather
than in one launcher's `.nextflow/`. This matters once runs start from
a web server as well as developer machines. Keep dev and prod apart,
as with the pools: e.g. `cache/wf5/test/` and `cache/wf5/prod/`, set
by the launcher alongside the pool. Pin the plugin version (rule 11).
Note in docs/azure that a cache entry is only useful while its work
dir still exists. The 14-day `workdata/` lifecycle rule therefore caps
how far back `-resume` can reach, and the `cache/` container has no
lifecycle rule of its own. Decide whether it needs one, and record
the decision.

### 3.4 Inputs and outputs

- **Reads** come in either of task 56's two modes. **Mode A**: a
  local `--data_dir` with relative paths. Nextflow uploads only the
  listed FASTQs, once per run, into the work dir's `stage-*` area
  under `az://workdata/wf5/`, and only `PARSE_SAMPLESHEET` receives
  them. **Mode B**: `az://` URIs in the sheet and no `--data_dir`, so
  nothing is uploaded. `--data_dir` may also be an `az://` prefix with
  relative paths, as the `--integration` preset uses.
- **`--samplesheet`** may be local. It is small and staged as normal.
- **`--outdir`** defaults to a local directory on the launcher, so the
  run's reports are directly inspectable. Every `publishDir` in this
  pipeline uses `mode: 'copy'`, which works from a blob work dir. An
  `az://` outdir must also work, and the script should accept one.

### 3.5 Credentials and config

One set of environment variable names, used by `.env.azure`,
`.env.sample`, `batch-helpers.sh`, `conf/azure.config`, the launcher and
docs/azure. Pick the names `batch-helpers.sh` and docs/azure already
use (`AZURE_BATCH_ACCOUNT_NAME`, `AZURE_BATCH_ACCESS_KEY`,
`AZURE_BATCH_ENDPOINT`, `STORAGE_ACCOUNT_STD`,
`AZURE_STORAGE_ACCOUNT_KEY`) and make `conf/azure.config` read those.
`conf/azure.config` must fail clearly when a required variable is unset.
A null `accountKey` silently produces an opaque auth error deep in the
plugin. Remove `.env.sample` entries nothing reads.

Pin both versions (rule 11): Nextflow `25.04.0`, exported as `NXF_VER`
by the launcher (26.x fails to parse `main.nf`, see the remote-run
skill), and the `nf-azure` plugin, pinned in the azure profile.

`conf/azure.config` must not set reference params. Those belong to the
launcher (§4.4), so one profile serves any bundle version.

## 4. Work items

1. **Pool templates**: production and test, per §3.1. Create
   `wf5_test` first and run every verification in §5 against it.
   Create `wf5` once §5 passes.
2. **Reduce the start task** to the `azcopy` install, per §3.2, and
   upload it to `scripts/wf5/setup.sh`. Set
   `copyToolInstallMode = 'off'` in `conf/azure.config`.
3. **Extracted refdata on blob**: §3.2. Publish `v2026.09_1`, the
   version CI currently pins, and verify completeness.
4. **`deploy/azure/run-wf5-azure.sh`**: replaces
   `deploy/azure/run-wf5.sh`, which is deleted. Behaviour (shell is
   fine; promote anything over trivial into a function):
   - Load `.env.azure` and verify required vars. Use `NXF_VER=25.04.0`.
   - Args: `--samplesheet`, `--data-dir` (optional — omit for a
     mode-B URI sheet; passed to Nextflow as **`--data_dir`**, since
     the kebab form silently falls through),
     `--outdir`, `--refs-version` (required, no hard-coded default that
     goes stale), `--pool` (default `wf5_test`; production must be asked
     for explicitly), `--polish`, `-resume`, plus
     two presets:
     - `--stub`: `-profile azure,stub -stub-run`, the cheap
       wiring/auth/pool check.
     - `--integration`: `tests/integration/samples.csv` against
       `az://test-data/wf5/<fixtures version>` (default matches
       `tests/integration/fetch_fixtures.sh`). Runs
       `tests/integration/assertions.sh <outdir>` on success, and the
       script's exit code reflects the assertions.
   - Derive all five reference params from `--refs-version`:
     `organelle_refs`, `protein_panel`, `blast_db`, `annotate_refs`,
     `refs_manifest`, all under `az://refdata-wf5/<version>/`.
   - Preflight: the pool exists and is `active`, and the `workdata`
     container is reachable. When `--data-dir` is local, print the
     total size of the files the sheet lists. Mode A uploads those over
     the launcher's link, so the cost should be visible up front.
     Fail on anything required that is missing (rule 18). Never fall
     back to a default bundle or pool.
   - Select the cache path from the pool (§3.3: test vs prod).
   - Work dir `az://workdata/wf5/` with a per-run subdirectory, so
     concurrent runs and post-run cleanup don't collide. Print it, along
     with the `az_fetch_nf_logs` command to pull task logs on failure.
   - Keep the existing interactive confirmation, with `--yes` to skip.
   - Print the path to the rendered `run-report.html` at the end.
5. **`batch-helpers.sh`**: update anything that hard-codes the old
   script location or SAS flow. Keep changes minimal.
6. **Docs**: rewrite docs/azure/README.md and docs/azure/setup.md for
   the new flow. Cover the `workdata/wf5` and `cache/wf5` prefixes,
   the extracted-bundle publish, the two pools and when to use each,
   the reduced start task and its new blob path, the token-renewal
   table (refdata token row goes, setup-script row stays), and the new
   launcher.
   Fix setup.md's dangling "§5.3" reference.
7. **Spec**:
   - spec/04-reference-data.md §4.4–4.5: the extracted bundle on blob is
     a published form alongside the tarball, with the same immutability
     and completeness rules.
   - spec/04a-azure-blob-storage.md: mark which prefixes are now used
     and by what.
8. **tasks/todo.md**: remove the Azure half of the `plant_mt.mmi` RSS
   item once §3.1's measurement is recorded. Leave the CI-runner half.

## 5. Verification

### 5.1 Local regression

`-profile stub,docker -stub-run` and `-profile integration` still pass.
This task should touch no pipeline code beyond config.

### 5.2 Azure stub run

`run-wf5-azure.sh --stub --refs-version v2026.09_1` completes on
`wf5_test`. It proves auth, pool scheduling, the start task's
`azcopy` being visible to task containers with Nextflow's install
turned off, container pulls on the node, `az://` reference staging
and publishing back. Fix wiring here before spending
on real tools.

### 5.3 Azure integration run

`run-wf5-azure.sh --integration --refs-version v2026.09_1` completes,
and `tests/integration/assertions.sh` passes against the published
output. That is the same bar as nightly CI, on different
infrastructure. Run with `-with-trace` and record in this task's
outcome notes:

- wall time per process and total
- reference staging time for `RECRUIT` and `VALIDATE`/BLAST tasks.
  This is the evidence behind §3.2's per-task staging choice, so state
  plainly whether it is material.
- peak RSS of the `plant_mt` processes in §3.1

Render and inspect `run-report.html` and one per-sample `report.html`.
Give the user their paths, per CLAUDE.md.

### 5.4 Read-source modes

Both of task 56's modes, run on the integration fixtures with
`assertions.sh` green each time:

- **Mode A, staged:** download the fixtures locally (as
  `tests/integration/fetch_fixtures.sh` does) and run with a local
  `--data-dir`. Confirm from the work dir listing that only the three
  listed FASTQs appear under `stage-*`, uploaded once. Confirm from the
  trace that raw reads are inputs to `PARSE_SAMPLESHEET` only.
- **Mode B, URIs:** a copy of `tests/integration/samples.csv` whose
  `reads` are full `az://test-data/wf5/<fixtures version>/…` URIs, run
  with no `--data-dir`. Confirm nothing is uploaded to `stage-*`.

The `--integration` preset (relative paths + `az://` `--data_dir`)
covers the third combination. Mode B can be a one-off manual check,
with the URI sheet kept under `tests/manual/`. It doesn't need its own
preset.

### 5.5 Negative checks

- Unset `AZURE_BATCH_ACCESS_KEY` → the launcher fails before Nextflow
  starts, naming the variable.
- `--refs-version` naming a version with no extracted bundle → the run
  aborts in `validateParams()` with the existing "bundle is missing or
  incomplete" message. Not a container-level failure.
- `--pool` naming a non-existent pool → preflight failure.
- `-resume` of the integration run from a **second** launcher checkout
  (or after deleting the local `.nextflow/`) reuses cached tasks via
  `cache/wf5/test/`.

### 5.6 Cleanup

Remove the run's `workdata/wf5/<run>/` prefix after a successful
verification run. Leave the pool scaled to zero.

## 6. Acceptance criteria

- [ ] Pools `wf5_test` and `wf5` created from their templates on a
      supported image, sized per §3.1.
- [ ] Start task installs a pinned `azcopy` only, from
      `scripts/wf5/setup.sh`, where task containers can see it. Refdata
      staging and its SAS placeholders are removed.
      `copyToolInstallMode = 'off'`.
- [ ] `conf/azure.config` targets the named pool, `workDir` under
      `az://workdata/wf5/`, cloud cache under `az://cache/wf5/`,
      consistent env var names, pinned plugins.
- [ ] Extracted `v2026.09_1` bundle published under
      `refdata-wf5/v2026.09_1/` and verified complete.
- [ ] `deploy/azure/run-wf5-azure.sh` exists. `deploy/azure/run-wf5.sh`
      is removed.
- [ ] Azure stub run passes.
- [ ] Azure integration run passes `tests/integration/assertions.sh`.
      Timings and RSS are recorded.
- [ ] Both read-source modes pass on Azure (§5.4). Mode A uploads only
      the listed reads.
- [ ] Negative checks in §5.5 behave as described.
- [ ] docs/azure and spec 04 / 04a updated. tasks/todo.md item trimmed.

## 7. Decisions recorded

- **Reference data:** staged per task from `az://`, not pre-staged on
  nodes (§3.2). Confirmed with the user 2026-10-12.
- **Start task:** kept, but only to install `azcopy`. Tool installation
  belongs to the pool, not the workflow, so Nextflow's install is
  turned off (§3.2). Confirmed with the user 2026-10-12.
- **Cache:** blob cloud cache under `cache/wf5/`, split test/prod
  (§3.3). Confirmed with the user 2026-10-12.
- **Test pool:** kept, so dev runs stay separate from prod (§3.1).
  Confirmed with the user 2026-10-12.
