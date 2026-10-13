# Task 67 — `ps` in the task images, so runs can be traced

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** open. No preconditions. Overlaps with task
61_intake_integrity_and_failure_visibility.md on two container lines
(§2.2). Whichever task lands second drops them from its diff.

## Overview

Nextflow's execution report (`-with-report`), timeline
(`-with-timeline`) and trace (`-with-trace`) all turn on per-task
metrics. Under Nextflow 25.04.0, the task wrapper's `nxf_trace_linux`
(`nextflow/executor/command-trace.txt` in the framework jar) then runs
this first:

```
command -v ps &>/dev/null || { >&2 echo "Command 'ps' required by nextflow to collect task metrics cannot be found"; exit 1; }
```

It runs **before** the task script. A process whose image has no `ps`
therefore fails outright, and does not just lose its metrics. With
`errorStrategy = 'ignore'`, that sample then drops out of the run
silently. No config option gives a timeline or report without metrics.

For this reason, task 57 made tracing opt-in in
`deploy/azure/run-wf5-azure.sh` (`--trace`). Neither the `azure-run`
nor the `remote-run` skill writes a report or timeline. So the
`tests/manual/output/` runs have no execution report, no timeline and
no per-process resource figures. That blocks:

- Reviewing a run's wall time and failures from one HTML page, which
  the user has asked for on every manual run.
- Sizing Batch pools and process labels from evidence.
- The open item from task 57: the `plant_mt` processes' peak RSS on
  Azure. `plant_mt.mmi` is 536 MB since task 28, and it is loaded by
  `RECRUIT`, `COVERAGE_GATE` and `BIN_TARGET`. Task 57 showed it fits
  on a 64 GiB node, but not by how much (tasks/todo.md, Reference
  data).

## 1. Survey (done 2026-10-14)

Every image in `conf/containers.config` was checked locally with
`docker run --rm --user $(id -u):$(id -g) --entrypoint "" <image> sh -c 'command -v ps'`:

| Image | Processes | `ps` |
|---|---|---|
| `neoformit/daff-wf5-scripts@sha256:2fcdfb…` (684d981) | `COLLATE`, `RUN_REPORT` | **missing** |
| `neoformit/daff-wf5-scripts@sha256:4fc343…` (5596079) | `ALLOCATE_GRAPH_SENTINELS`, `ANNOTATE_GRAPH_SVG`, `SELECT_GENETIC_CODE`, `EXTRACT_BARCODES`, `ORGANELLE_MAP`, `TRANSLATE_ANNOTATION_CDS`, `JOIN_ANNOTATION_SCORES` | **missing** |
| `python:3.12-slim` (by tag, not digest) | `VALIDATE_SAMPLESHEET`, `PARSE_SAMPLESHEET`, `MEDAKA` (placeholder) | **missing** |
| `neoformit/daff-wf5-coverage-gate` | `COVERAGE_GATE`, `RECRUIT` | `/usr/bin/ps` |
| `neoformit/daff-wf5-bin-target` | `BIN_TARGET` | `/usr/bin/ps` |
| biocontainers: nanoplot, chopper, filtlong, flye, bandage_ng, blast, mitos, miniprot | tool stages | `/usr/bin/ps` |

So the fix is narrower than first thought. **No biocontainer and no
mulled image needs to change.** Only `scripts/Dockerfile` lacks
`procps`, plus the remaining `python:3.12-slim` placeholders. The
mulled images got `procps` through their conda base, and the
biocontainers have a real (procps) `ps`, not BusyBox's.

## 2. Scope

### 2.1 `scripts/Dockerfile`: add `procps`

Install it in the runtime layer, before the `pip install`:

```dockerfile
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
```

Add `ps --version` to the image's existing smoke-test `RUN`, so a build
without it fails. The `TEST=1` layer purges only `gcc zlib1g-dev`, so
it does not remove `procps`. Check `apt-get autoremove` does not take
it either.

`images.yml` rebuilds on a push to `main` that touches the Dockerfile.
Resolve the new digest as in spec/containers.md §1 and pin it.

### 2.2 `conf/containers.config`: one scripts image, and no `python:3.12-slim`

- Move every `neoformit/daff-wf5-scripts` entry, both the `684d981` and
  the `5596079` digests, to the new digest. Two pinned versions of one
  image is drift with no reason behind it.
- Replace `python:3.12-slim` for `VALIDATE_SAMPLESHEET` and
  `PARSE_SAMPLESHEET` with the new scripts digest. Task 61 §1 plans the
  same change for rule 12 (digest pinning). If task 61 has landed
  first, these lines are already done.
- `MEDAKA` is a `touch` placeholder (task 17 is shelved). Move it to the
  scripts digest too, so that `--polish` runs can be traced. Keep its
  `TODO P2` comment naming the real medaka biocontainer.
- Remove the P0 header comment about `python:3.12-slim` and the
  `TODO: replace python:3.12-slim` comment once no entry uses it.

### 2.3 Report, timeline and trace on by default

- **`nextflow.config`:** set `report.overwrite`, `timeline.overwrite`
  and `trace.overwrite` to `true`. Without that, a `-resume` into the
  same outdir aborts because the file already exists. Setting these does
  not enable the reports; the flags or the launcher do that.
- **`deploy/azure/run-wf5-azure.sh`:** always pass
  `-with-report "${OUTDIR}/execution_report.html"`,
  `-with-timeline "${OUTDIR}/execution_timeline.html"` and
  `-with-trace "${OUTDIR}/execution_trace.txt"`, for every mode,
  `--integration` included. Remove the `--trace` flag and its comment,
  and update the usage header and docs/azure/README.md, which still
  says `--trace` needs `ps`. The file names are the defaults. Pick them
  once and use them in both skills (§2.4).
- Confirm that `tests/integration/assertions.sh` ignores the three new
  files in the outdir. It checks per-sample bundles, so it should, but
  check for any whole-outdir listing.

### 2.4 Final step: skills

Do this last, after an Azure run and a `daff-admin` run have both
produced the three files (§3).

- **`.claude/skills/azure-run/SKILL.md`:** the launcher now writes the
  files, so no extra flags are needed. In "Review", point the user at
  `execution_report.html` and `execution_timeline.html` next to
  `run-report.html`. Add to "If a task fails" that the report's tasks
  table (exit status, work dir) is the quickest route to the failing
  task, before `az_fetch_nf_logs`.
- **`.claude/skills/remote-run/SKILL.md`:** add these flags to the
  step 4 `nohup` command:
  ```
  -with-report   /mnt/data/wf5/output_<run-name>/execution_report.html \
  -with-timeline /mnt/data/wf5/output_<run-name>/execution_timeline.html \
  -with-trace    /mnt/data/wf5/output_<run-name>/execution_trace.txt \
  ```
  Step 6's pull of the outdir brings them back with everything else.
  Mention them in step 6's review the same way.
- In both skills, say that a run whose task fails with "Command 'ps'
  required by nextflow…" means an image has lost `procps`, and that
  this is a regression of this task, not a sample failure.

## 3. Testing

- **Image:** `docker run --rm --user $(id -u):$(id -g) <new digest> ps --version`
  succeeds, and `scripts/pytest.sh` passes in the `TEST=1` image.
- **Local stub:** `-profile stub,docker -stub-run -with-report … -with-timeline … -with-trace …`
  completes with no ignored task. Before this task, every
  scripts-image process fails this way, so that is the regression
  check.
- **Local integration:** `-profile integration,docker` with the three
  flags. `assertions.sh` is green, and `execution_trace.txt` has a
  non-empty `peak_rss` for every process.
- **Azure:** `deploy/azure/run-wf5-azure.sh --integration --refs-version v2026.09_1`
  writes the three files to the outdir, `assertions.sh` passes, and no
  task is ignored. Read the `plant_mt` (`INT-PLANT-01-mt`) peak RSS for
  `RECRUIT`, `COVERAGE_GATE` and `BIN_TARGET` from the trace. Close the
  task 57 item in tasks/todo.md and the Azure half of the
  `plant_mt.mmi` item.
- **`daff-admin`:** one real sample through the updated `remote-run`
  step 4 command returns the three files.

## 4. Spec updates

- spec/01-pipeline-flow.md §1a: every task image must provide a procps
  `ps` (not BusyBox), because metrics collection runs it before the
  task script and fails the task without it.
- spec/containers.md §1: note `procps` in the scripts image. Update the
  "As of task 30" paragraph, which still lists `python:3.12-slim`
  placeholders.

## 5. Acceptance criteria

- [ ] `scripts/Dockerfile` installs `procps`, and the smoke test checks it.
- [ ] `conf/containers.config` has one `daff-wf5-scripts` digest and no
      `python:3.12-slim`.
- [ ] The launcher writes the report, timeline and trace for every mode,
      and `--trace` is gone.
- [ ] A `-resume` into an existing outdir does not abort on an existing
      report file.
- [ ] Local stub, local integration and Azure integration all pass with
      metrics on, with no ignored task.
- [ ] The `plant_mt` peak RSS is recorded and the tasks/todo.md items
      are closed.
- [ ] Both skills are updated (§2.4), after the runs above.
