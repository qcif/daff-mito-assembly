# Task 67 — `ps` in the task images, so runs can be traced (stub)

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** stub. Scope it fully before starting. No preconditions.
Best done after task 61_intake_integrity_and_failure_visibility.md,
which moves the samplesheet processes off `python:3.12-slim`, leaving
fewer images to fix.

## Overview

Nextflow's trace report (`-with-trace`) records each task's wall time,
CPU and peak memory (RSS). It collects these by running `ps` inside the
task's container. Most of the images this pipeline uses do not contain
`ps`. Run with tracing on, those tasks lose their metrics or fail. So
task 57 made tracing opt-in in `deploy/azure/run-wf5-azure.sh`
(`--trace`) and left it off by default.

Without the trace there are no per-process resource figures from Azure.
That blocks two things:

- Sizing Batch pools and process labels from evidence.
- The open item from task 57: the `plant_mt` processes' peak RSS on
  Azure. `plant_mt.mmi` is 536 MB since task 28, and it is loaded by
  `RECRUIT`, `COVERAGE_GATE` and `BIN_TARGET`. Task 57 showed it fits
  on a 64 GiB node, but not by how much (tasks/todo.md, Reference
  data).

## Known so far

- The bespoke images come from `scripts/Dockerfile` (`python:3.12-slim`,
  no `procps`) and the mulled images `daff-wf5-coverage-gate` and
  `daff-wf5-bin-target`. All are listed in `conf/containers.config`.
- The biocontainers images (NanoPlot, chopper, filtlong, Flye,
  BandageNG, BLAST, MITOS2, miniprot, and barrnap from task
  65_plastid_noncds_barrnap.md) are built on a minimal base. Whether
  each has a `ps` that Nextflow's trace accepts (BusyBox `ps` may not)
  has not been checked.

## Scope to decide when this is picked up

1. **Survey.** For each image in `conf/containers.config`, check for
   `ps` and whether a short traced local run of that process records
   `peak_rss`. Make a table.
2. **Bespoke images:** add `procps` to `scripts/Dockerfile`. This
   triggers `images.yml`, so update the digests in
   `conf/containers.config`. For the two mulled images, check whether
   `mulled-build` can add `procps` from conda-forge, rather than
   hand-writing a Dockerfile (spec §1a container order).
3. **Biocontainers without a usable `ps`:** pick one approach for all
   of them, and record the provenance. The options are a mulled
   rebuild with `procps` added, or leaving them untraced and getting
   the figures another way, for example Azure Batch node metrics. Do
   not hand-patch a pinned upstream image.
4. **Launcher:** make `--trace` the default once every image is traced.
5. **Measure** `plant_mt` peak RSS on Azure with the `INT-PLANT-01-mt`
   fixture. Close the task 57 item in tasks/todo.md, and the Azure half
   of the `plant_mt.mmi` item.

## Testing (outline)

- Each changed image passes its existing tests and gives the same
  outputs. Run the integration profile locally and confirm the
  assertions are green with digests bumped.
- A local `-with-trace` integration run has non-empty `peak_rss` for
  every process.
- An Azure stub or integration run with `--trace` gives the same.

## Spec updates (outline)

- spec/01-pipeline-flow.md §1a: task images must provide `ps`
  (procps) for tracing.
- spec/containers.md: a note per image.
