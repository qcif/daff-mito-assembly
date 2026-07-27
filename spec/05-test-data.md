## 5. Test data

Per kingdom, one clean + one contaminated dataset. Skim depth
(~2× nuclear-equivalent) throughout — see [brief.md §2](brief.md).

| Kingdom | Clean target | Minor-contamination scenario |
|---------|-------------|------------------------------|
| plant | known reference plant (e.g. *Arabidopsis*) ONT skim | plant + trace insect contaminant |
| animal | known insect (e.g. *Drosophila*) ONT skim | insect + trace plant host DNA |

Acceptance criteria per dataset:
- Clean: full organelle assembled, all panel loci extractable with valid ORF,
  full annotation produced by `ANNOTATE`.
- Contamination: target dominant assembly selected; any low-coverage secondary
  contigs recorded in diagnostics; target assembly unaffected by trace
  non-target reads.

## 5a. Tests

Two testing surfaces, one per language, both run in CI on every push
inside the same containers the workflow uses at runtime
([§1a](01-pipeline-flow.md#1a-engineering-constraints)):

- **Python (`bin/*.py`, `wf-report-boilerplate/`) — `unittest`.**
  One test file per module under `tests/unit/`, imported directly (no
  Nextflow harness). Each [custom-logic component](02-stages.md#22-custom-logic-components)
  (C1–C7) gets its own test module. Fixtures live under
  `tests/fixtures/` — synthetic GFA files, minimal FASTA snippets,
  hand-authored `seqkit stats` outputs — small enough to check in.
  `flake8` runs alongside the tests (per user CLAUDE.md).
- **Nextflow (`main.nf`, `modules/local/*.nf`) — [nf-test](https://www.nf-test.com/).**
  Per-module tests exercise each process's channel wiring, `input:` /
  `output:` shapes, and stub behaviour in isolation. A workflow-level
  test replays the P0 `-profile test` end-to-end run and asserts the
  per-sample output layout in [§0](00-overview.md#0-input-sample-sheet). Reference
  bundles are pointed at the same `tests/fixtures/` used by the Python
  unit tests.

**Coverage target: 100%.** Where a component is genuinely untestable in
CI — because it is variable (calls out to a network service, is
non-deterministic), expensive (needs a real ONT dataset or a full
BLAST DB), or too complicated to fixture (`METAFLYE`, `MEDAKA`,
`ANNOTATE`, `ORGANELLE_MAP`) — mock the tool call at the process boundary and
test only our wrapper logic. Record every mocked-out surface in the
per-component test module's docstring so it's visible what is *not*
being covered.

**Rule of thumb for what to mock:**

- Off-the-shelf tools invoked by pure-wrapper processes
  (`NANOPLOT_*`, `CHOPPER`, `FILTLONG`, `METAFLYE`, `MEDAKA`,
  `BANDAGE_NG`, `BLAST_VALIDATE`, `ANNOTATE`, `ORGANELLE_MAP`,
  `MINIPROT_EXTRACT`) — mock at nf-test level; assert we assemble the
  command line correctly and consume the output shape correctly. Do
  not try to actually run the tool.
- Custom-logic components (C1–C7 in [§2.2](02-stages.md#22-custom-logic-components))
  — real Python unit tests with real fixtures. These are our code and
  they must be covered.
- Network fetches (`scripts/build_refs.sh` calls to NCBI FTP) — mock
  the fetch, real code for parsing.

CI matrix reports coverage per component; a drop below 100% for
non-mocked lines fails the build.

## 5b. CI

CI runs on **GitHub Actions**. Two concerns: run every test surface
defined in [§5a](#5a-tests), and publish the bespoke container images
listed in [§2.2](02-stages.md#22-custom-logic-components) so the workflow has
something to `docker pull`.

**Workflows** (one YAML file each, under `.github/workflows/`):

- **`tests.yml`** — runs on every push and PR to any branch.
  1. `nf-core lint` (or `nextflow lint main.nf` if the nf-core scaffold
     is skipped — see [tasks/1_scaffold.md §9 q1](../tasks/1_scaffold.md)).
  2. **Container-coverage check** — assert every process in
     `modules/local/` resolves to a `container` directive via
     `conf/containers.config`, and that no tag is `latest`
     ([§1a](01-pipeline-flow.md#1a-engineering-constraints)). A process without a
     container fails CI.
  3. **Stageable-file check** — grep-based lint asserts no process
     `script:` block contains a bare `${params.<file>}` interpolation
     that isn't wrapped in `file(...)` or on the scalar-param allow-list
     ([§1a](01-pipeline-flow.md#1a-engineering-constraints)).
  4. **Python** — `unittest` on `bin/*.py` + `wf-report-boilerplate/`
     with `flake8`; each component's tests run inside its own
     container image ([§2.2](02-stages.md#22-custom-logic-components)) so the
     runtime environment matches production exactly.
  5. **Nextflow** — `nf-test` per-module + workflow-level; the
     workflow-level test is the P0 `-profile test -stub` end-to-end
     pass. Container runtime enabled so image-pull failures surface
     here.
  6. **Coverage report** — per-component coverage published as a job
     summary; <100% on non-mocked lines fails the build.

- **`images.yml`** — rebuilds bespoke container images. Triggered
  **only** on push to `main` (or on a release tag) touching one of:
  - `scripts/requirements.txt` (or per-image `containers/<name>/requirements.txt`)
  - `scripts/Dockerfile` (or per-image `containers/<name>/Dockerfile`)
  - the workflow file itself.

  **Not triggered by Python source changes** — images bundle the
  interpreter + pinned dependencies only; the source under
  `scripts/` / `bin/` is delivered to the container at runtime (see
  [§2.2](02-stages.md#22-custom-logic-components) container pattern). Editing a
  Python file rebuilds nothing, ships nothing, and iterates in
  seconds.

  For each of the five bespoke images from
  [§2.2](02-stages.md#22-custom-logic-components) — `wf5/samplesheet`,
  `wf5/coverage-gate`, `wf5/bin-target`, `wf5/barcode-validate`,
  `wf5/report`:
  1. Build with `docker/build-push-action`.
  2. Tag with the short git SHA of the commit that touched deps + (on
     release) the semver tag.
  3. Push to **DockerHub** — free for public repos, need to provide DockerHub
     API credential for image push.
  4. Emit a build-provenance attestation (`actions/attest-build-provenance`)
     — biosecurity reporting benefits from a verifiable chain from
     source to image.
  Path filters are per-image so a change to only `wf5/coverage-gate`'s
  requirements doesn't rebuild the other four.

- **`ref-build.yml`** *(optional, added when the reference-build script
  in [§4.4](04-reference-data.md#44-consolidated-build-script) is wired)* — scheduled
  monthly, invokes `scripts/build_refs.sh` in the `wf5/ref-build` image,
  publishes the resulting `refs/v<YYYY.MM>/` bundle to a release
  artefact (or object store). Manual `workflow_dispatch` trigger for
  ad-hoc rebuilds.

**Branch protection:** `main` requires `tests.yml` green before merge.
`images.yml` runs post-merge so the image tag matches the merged
commit.

**Local re-run:** every check in `tests.yml` must be runnable outside
CI (`nf-test`, `python -m unittest`, `flake8`, and the two grep
checks). No CI-only tooling.
