# Task 56 — Samplesheet read sources: staged data dir or blob URIs

**Phase:** P5 prerequisite — remote-executor portability (CONSTITUTION
rule 17).
**Goal:** Support two ways of supplying reads, and make both work on
Azure Batch without copying input data into more processes than
consume it:

- **Mode A — CLI / development.** The sheet lists relative paths and
  `--data_dir` names a local directory, as today. Nextflow uploads only
  the listed FASTQs into the blob work dir. Only the process that reads
  raw FASTQ (`PARSE_SAMPLESHEET`) receives them.
- **Mode B — web server.** The user uploads reads directly to Azure
  Blob. The service writes a sheet whose `reads` entries are `az://`
  URIs, and `--data_dir` is not passed.

**Blocks:** task 57_azure_batch_deployment.md.

**Related spec:** spec/00-overview.md §0 (samplesheet schema and
validation rules), spec/01-pipeline-flow.md §1a (staging),
spec/02-stages.md §2.2 (C1). brief.md §4 (Inputs) also states the
"absolute paths rejected" rule. It sits under the brief's *flexible*
constraints and must be revised to match.

---

## 1. Overview

Every run starts from a samplesheet. Each row lists the sample's
nanopore FASTQ files. Today those must be paths relative to a data
directory (`--data_dir`). The pipeline handles them like this:

1. `VALIDATE_SAMPLESHEET` takes the CSV **and the whole data directory**
   as inputs.
2. Inside its container, `bin/parse_samplesheet.py` joins each relative
   path onto the data directory, calls `Path.resolve()`, checks the file
   exists, and writes the resulting **absolute path** into
   `samples.normalised.json`.
3. `main.nf` turns each absolute string back into a file with `file(it)`
   and hands it to `PARSE_SAMPLESHEET`, which concatenates a sample's
   files into one FASTQ. Every later step (`NANOPLOT_RAW`, `CHOPPER`,
   …) consumes that concatenated output, never the raw files.

On one Docker host this works by accident. On Azure Batch it breaks.
Batch nodes share no filesystem, so:

- the `path data_dir` input copies the **entire** data directory onto
  the node running one small validation task;
- the resolved paths are node-local
  (`/mnt/batch/tasks/.../wd/data/x.fastq.gz`) and vanish with the task;
- `file()` on the launcher then points at nothing, and reads cannot be
  staged into `PARSE_SAMPLESHEET`.

There is also no way to name reads that already live in blob storage,
which is how the planned web front end will supply them.

The fix is to resolve each `reads` entry **on the launcher**, in
`main.nf`, into a real Nextflow `Path`. That path is either
`file(params.data_dir).resolve(rel)` or `file('az://…')`. Nextflow then
moves each file only where it is declared as an input (§2.3). The
Python validator keeps every rule it can check from the CSV alone.

## 2. Design

### 2.1 What a `reads` entry may be

Each pipe-delimited entry is classified as one of:

| Kind | Example | Accepted when |
|---|---|---|
| relative path | `run17/A.fastq.gz` | `--data_dir` is set |
| blob URI | `az://reads-upload/job-123/A.fastq.gz` | always |
| anything else | `/data/A.fastq.gz`, `~/A.fq`, `https://…`, `s3://…` | never: rejected with an actionable error |

Rules, all enforced by C1 (`parse_samplesheet.py`) so they are
unit-tested and reported together with every other sheet error:

- **Local absolute paths stay rejected.** They were rejected so the
  sheet doesn't depend on the machine that wrote it, and that reason
  still holds. A URI is the one exception, because it names a location
  every executor can reach.
- **Allowed URI schemes are an explicit allow-list**, `az://` only for
  now, declared as a constant in `parse_samplesheet.py`. Adding
  `s3://` or `gs://` later is a one-line change plus a test. Do not
  accept "any scheme": a typo'd or unsupported scheme must fail at
  parse, not deep inside a staging error (rule 18).
- **A sheet is all-relative or all-URI, never mixed.** Each mode
  produces one kind, and no current use case needs both. A mixed sheet
  most likely means a hand-edit gone wrong. Reject it with an error
  naming the mode each kind implies. This is a recommendation; if the
  implementer finds a real use for mixing, drop the rule and record why.
- **Relative entries without `--data_dir` fail** at C1, listing the
  offending rows.
- **`--data_dir` with an all-URI sheet** produces a warning, not an
  error: the flag is harmless but unused. It is not silently ignored.
- The existing extension check applies to both kinds (judged on the
  path part of a URI). So does the empty-entry check.

C1 needs to know whether `--data_dir` was supplied, not its contents.
Pass it a boolean flag (e.g. `--data-dir-set`) rather than the
directory. Nothing in C1 should touch the filesystem for reads any
more.

`samples.normalised.json` keeps its shape. `reads` holds each entry
**verbatim**: relative paths stay relative, URIs stay URIs. Optionally
it also carries a sheet-level `read_source: relative|uri` field for
`main.nf` and the run manifest. Record the contract change under C1 in
spec §2.2.

### 2.2 Resolution and existence check move to the launcher

`main.nf` turns each entry into a `Path`. Relative entries become
`file(params.data_dir).resolve(rel)` and URIs become `file(uri)`. It
checks existence **once, holistically**: one error listing every
missing file, before any per-sample work starts. That matches C1's
current behaviour and spec §0's "missing files fail the whole run — no
partial execution".

Pseudocode (illustrative only):

```
rows = VALIDATE_SAMPLESHEET.out.json | splitJson | toList
rows.map { all ->
    missing = []
    for row in all, entry in row.reads:
        f = toPath(entry)          // data_dir.resolve(entry) or file(entry)
        if not f.exists(): missing << "${row.sample_id}: ${entry}"
    if missing: error "reads not found:\n" + missing.join("\n")
    return all
}
| flatten-to-rows
| map { row -> [ meta(row), row.reads.collect { toPath(it) } ] }
```

Choose whatever channel shape is clearest. Keep the per-sample tuple
identical in shape to today's.

`az://` existence checks run on the launcher with the Azure
credentials from the azure profile. Without them, e.g. a URI sheet run
under `-profile docker` with no Azure config, `exists()` fails. The
error must say that credentials are missing, not that the reads are.
Check what Nextflow 25.04 actually raises and wrap it if it is opaque.

`validateParams()`:

- `--data_dir` becomes optional. Its existence check runs only when it
  is set.
- Whether it is *required* is C1's call (§2.1), because only C1 knows
  whether the sheet holds relative entries. Do not duplicate that rule
  in Groovy.
- `--data_dir` may itself be an `az://` prefix, i.e. relative entries
  resolved against a blob prefix. This falls out of using `file()` and
  costs nothing. Task 57's `--integration` preset relies on it. Confirm
  `file('az://container/prefix').exists()` is true for a virtual
  directory with no blob at the prefix itself. If it isn't, record that
  and do the check another way.

### 2.3 Staging: who receives the reads, and when

There is no pipeline code to write here. This section records the
behaviour the design relies on, and §4 verifies it.

- **`VALIDATE_SAMPLESHEET` takes no `data_dir` input.** This one change
  removes the "copy the whole directory to a node" cost.
- **Only listed files move.** Each read becomes its own `Path`, so
  Nextflow never sees the directory. Unlisted files in `--data_dir` are
  never uploaded.
- **Mode A on a blob work dir:** a local read is a *foreign* file to
  the executor. Nextflow uploads it once per run into
  `<workDir>/stage-<session>/…` and reuses that copy for every task
  that consumes it. Under task 57 that location is
  `az://workdata/wf5/…`, so the 14-day lifecycle rule on `workdata/`
  covers staged inputs too. Note the upload time is paid on the
  launcher's uplink, which is why mode A is the development mode.
- **Mode B:** URIs are already where the executor can read them, so
  nothing is uploaded.
- **Only `PARSE_SAMPLESHEET` consumes raw reads**, one task per sample,
  receiving only that sample's files. Everything downstream consumes its
  concatenated output, which lives in the work dir like any
  intermediate. Keep it that way. Add a comment at the
  `PARSE_SAMPLESHEET` call in `main.nf` saying it must remain the single
  raw-read consumer, and why.
- **Local Docker runs are unchanged.** A local work dir means no upload
  and Nextflow symlinks as before.

### 2.4 Provenance

`run_manifest.json` already snapshots the samplesheet, so URIs and
relative paths are recorded verbatim. Add the read-source mode (and
`--data_dir` when set) to the run manifest if it is cheap. An auditor
should be able to see where the reads came from without decoding the
sheet (rule 16). If this touches C7's schema
(`assets/run_manifest.schema.json`), bump the schema and test it like
any other C7 change.

## 3. Unit tests (`scripts/tests/`)

C1 targets 100% branch coverage (rule 14). Mock at the filesystem /
tool boundary per task 27. No test may depend on a real data directory
or a real blob.

- Relative entries are emitted verbatim (single and pipe-delimited).
- `az://` entries are emitted verbatim. Extension is checked on the
  URI's path.
- Rejected: local absolute paths, `~`, non-allow-listed schemes
  (`https://`, `s3://`), malformed `az://` (no container / no blob
  path).
- Mixed relative + URI sheet → one error naming both kinds.
- Relative entries with the flag unset → error listing the rows.
- URI-only sheet with the flag set → warning, success.
- Remove the tests that asserted resolved absolute paths and
  Python-side "file not found".
- If §2.4 lands, update C7 tests for the new manifest field.

Run `scripts/pytest.sh` and confirm 100% branch coverage on
`parse_samplesheet.py`. Run flake8 on changed files.

## 4. Integration / stub verification

No new fixtures. Expected outputs are unchanged, because the
integration profile is mode A with a local `--data_dir`.

- `-profile stub,docker -stub-run` passes. The stub branch of
  `VALIDATE_SAMPLESHEET` emits `"reads": []`, so the existence check
  iterates over nothing.
- `-profile integration` passes and `tests/integration/assertions.sh`
  is green.
- `.nextflow.log` / trace confirm that `VALIDATE_SAMPLESHEET` stages no
  data directory and that raw reads are inputs to `PARSE_SAMPLESHEET`
  only.
- Manual negative check: two bogus relative entries → one fatal error
  naming both, before any per-sample process starts.
- Mode B and staged mode A are verified on Azure in task 57 (§5.3 there).
  This task can't exercise `az://` without that environment.

## 5. Docs / spec

- spec/00-overview.md §0: the `reads` column accepts relative paths
  (with `--data_dir`) or `az://` URIs (without). Add the
  mixed-sheet rule, the scheme allow-list and the now-optional
  `--data_dir`. Rewrite the "absolute path not allowed" rule to keep
  rejecting local absolute paths while naming URIs as the alternative.
  Note there that the flag is `--data_dir` on the Nextflow command line.
- brief.md §4 Inputs: same revision, flagged in the commit as a
  flexible-constraint change.
- spec/02-stages.md §2.2 C1: reads stay verbatim, there is no
  filesystem access, and existence is checked on the launcher.
- spec/01-pipeline-flow.md §1a: samplesheet reads as a worked example
  of "resolve on the launcher, stage as `Path`", plus the
  single-raw-read-consumer rule from §2.3.

## 6. Acceptance criteria

- [ ] `VALIDATE_SAMPLESHEET` takes no `data_dir` input. C1 does no
      filesystem access for reads.
- [ ] Sheets of relative paths (with `--data_dir`) and of `az://` URIs
      (without) both validate. Mixed sheets, local absolute paths and
      non-allow-listed schemes are rejected with actionable errors.
- [ ] Missing reads fail the run once, before fan-out, listing every
      missing entry.
- [ ] `PARSE_SAMPLESHEET` is the only process receiving raw reads.
      This is confirmed from the trace and commented in `main.nf`.
- [ ] pytest green with C1 at 100% branch coverage. flake8 clean.
- [ ] Stub and integration runs pass unchanged.
- [ ] spec §0, §1a, §2.2 and brief.md §4 updated.
