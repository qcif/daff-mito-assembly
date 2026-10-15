# Container images — where to find them

Working notes for locating/rebuilding the images `conf/containers.config`
pins by SHA. Not a spec chapter with its own numbering; a lookup table,
kept up to date as images change.

## 1. Custom-logic image (`neoformit/daff-wf5-scripts`)

One shared image for every `bin/*.py` custom-logic stage (C1–C8) — deps in
the image, code staged at runtime via Nextflow's `bin/` auto-PATH
mechanism (spec §1a). Built from [`scripts/Dockerfile`](../scripts/Dockerfile)
+ [`scripts/requirements.txt`](../scripts/requirements.txt).

- **CI build:** [`.github/workflows/images.yml`](../.github/workflows/images.yml)
  — triggers only when `scripts/requirements.txt`, `scripts/Dockerfile`, or
  the workflow file itself changes (never on `bin/*.py` edits). Pushes to
  Docker Hub `neoformit/daff-wf5-scripts`, tagged `type=sha` (short commit
  SHA) + semver on release tags.
- **Manual build/push** (when a stage needs the image before CI has run,
  or for local testing):
  ```sh
  docker build -t neoformit/daff-wf5-scripts:<tag> scripts/
  docker push neoformit/daff-wf5-scripts:<tag>
  docker inspect --format='{{index .RepoDigests 0}}' \
      neoformit/daff-wf5-scripts:<tag>
  ```
- **Resolving an already-pushed tag's digest** (no rebuild needed):
  ```sh
  docker pull neoformit/daff-wf5-scripts:<tag>
  docker inspect --format='{{index .RepoDigests 0}}' \
      neoformit/daff-wf5-scripts:<tag>
  ```
- **As of task 67:** every process that uses the scripts image
  (including the `MEDAKA` placeholder) is pinned to one digest,
  `sha256:ca6b25e55d4459b24ab2d2140a377b81999f02bd9a874776f36a564aa6d7389b` (tag `9481373`). No process uses `python:3.12-slim`.
  The image installs `procps`, because Nextflow runs `ps` before the task
  script when collecting metrics (spec/01-pipeline-flow.md §1a); the
  Dockerfile smoke test runs `ps --version`.

## 2. Off-the-shelf biocontainers

Single-tool images from `quay.io/biocontainers` — no build step, just
resolve and pin the digest:

```sh
docker pull quay.io/biocontainers/<tool>:<tag>
docker inspect --format='{{index .RepoDigests 0}}' \
    quay.io/biocontainers/<tool>:<tag>
```

Search [biocontainers.pro](https://biocontainers.pro) or
`quay.io/biocontainers/<tool>` for available tags. See
[spec §1a](01-pipeline-flow.md#1a-engineering-constraints) for the
container-selection order (single-tool biocontainer → `mulled-build`
bundle → hand-written Dockerfile, in that preference order) and
[task 21 §1](../tasks/completed/21_blast_validate.md#1-container-pin--confcontainersconfig)
for a worked example.

## 3. `mulled-build` bundles

Used where a process needs more than one bioconda-packaged tool together
(e.g. `RECRUIT` / `COVERAGE_GATE` share `minimap2` + `seqtk`, mirrored to
`neoformit/daff-wf5-coverage-gate`). Composed with `galaxy-tool-util`'s
`mulled-build`, not a hand-written Dockerfile — see
[spec §1a](01-pipeline-flow.md#1a-engineering-constraints) for the
invocation shape. Record the exact `mulled-build` command used in the
task that introduces the bundle.

## 4. Checking what's already pinned

`conf/containers.config` is the single source of truth for what's
actually wired in; this file is a "how do I get one" reference, not a
duplicate inventory. `grep withName conf/containers.config` to see
current pins.
