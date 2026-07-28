# Task 12 — Restore fast CI after the stub-profile migration

**Phase:** P1 CI maintenance.  
**Goal:** Fix both deterministic failures in the `Tests` workflow introduced
or exposed by the move to stub-only CI, without changing pipeline behaviour.

**Failure investigated:** [GitHub Actions run 30327661131](https://github.com/qcif/daff-mito-assembly/actions/runs/30327661131),
for commit `a11a657d6efde18defd7a6bd9d28d78767ab855d`
(`Reduce CI nextflow test to -stub-run only`), on 2026-07-28.

## Diagnosis

The workflow has two independent failures.

### 1. `Python unit tests + flake8`

The `flake8` step still passes the deleted directory
`wf-report-boilerplate/`:

```text
wf-report-boilerplate/:0:1: E902 FileNotFoundError:
[Errno 2] No such file or directory: 'wf-report-boilerplate/'
```

Commit `50c60c3` moved that tree to
`scripts/report/boilerplate/`, but `.github/workflows/tests.yml` was not
updated. The current exclusion is stale for the same reason.

### 2. `Nextflow stub end-to-end`

The workflow requests both a Nextflow execution report and trace:

```text
-with-report tests/output/nextflow_report.html
-with-trace tests/output/nextflow_trace.txt
```

Those observers enable per-task metric collection. The first stub runs in
`python:3.12-slim`, which does not contain `ps`, so Nextflow fails
`VALIDATE_SAMPLESHEET` before the stub command can complete:

```text
Command 'ps' required by nextflow to collect task metrics cannot be found
```

This is an instrumentation/container compatibility failure, not a
samplesheet or channel-wiring failure. Resource-utilisation metrics from
`touch`/`echo` stub commands have no diagnostic value in fast CI.

The proposed change was checked against an isolated archive of the failing
commit: the same `-profile stub -stub-run` invocation, with only
`-with-report` and `-with-trace` omitted, completed all 17 processes
successfully in 11 seconds using the existing container configuration.

## Implementation

Make the following changes in `.github/workflows/tests.yml`.

### 1. Correct the lint inputs

Replace:

```yaml
- name: flake8
  run: flake8 scripts/ bin/ wf-report-boilerplate/ --max-line-length=100 --exclude=wf-report-boilerplate/static/
```

with:

```yaml
- name: flake8
  run: flake8 scripts/ bin/ --max-line-length=100 --exclude=scripts/report/boilerplate/static/
```

`scripts/` already recursively includes `scripts/report/boilerplate/`, so
do not add that directory a second time.

### 2. Keep the stub run free of resource profiling

Change the invocation to:

```yaml
- name: Run stub pipeline
  run: |
    nextflow run . -profile stub -stub-run \
      --outdir tests/output
```

Remove the `Upload Nextflow reports` step because the report and trace files
will no longer be produced.

Add a failure-only diagnostics upload instead:

```yaml
- name: Upload Nextflow failure diagnostics
  if: failure()
  uses: actions/upload-artifact@v4
  with:
    name: nextflow-failure-diagnostics
    path: |
      .nextflow.log
      work/**/.command.{sh,run,out,err,log,exitcode}
    if-no-files-found: ignore
    retention-days: 7
```

This preserves actionable logs when the stub workflow itself fails without
requiring metric collection inside every process container. Keep the
existing output-layout verification unchanged.

## Exit criteria

- The `Python unit tests + flake8` job passes on a clean checkout.
- The lint command covers Python under both `scripts/` and `bin/` and no
  longer references `wf-report-boilerplate/`.
- The `Nextflow stub end-to-end` job completes all 17 stub processes with
  Nextflow `25.04.0`.
- The existing output-layout check passes for `STUB-01`,
  `run_manifest.json`, and `run-report.html`.
- A successful run does not upload the failure-diagnostics artifact.
- A deliberately failing local test of one stub demonstrates that the
  failure-diagnostics step finds `.nextflow.log` and available
  `.command.*` files.
- The complete `Tests` workflow is green for both `push` and
  `pull_request`.

## Not in scope

- Installing `procps` into every current or future process container.
  Real-tool images should include the runtime requirements of their tools;
  stub CI should not impose a cross-container metrics dependency.
- Building or retagging `neoformit/daff-wf5-scripts`.
- Changing any Nextflow module, channel, stub body, output layout, or
  pipeline parameter.
- Implementing the nightly real-tool integration workflow from
  [task 11](completed/11_integration_tests.md).
