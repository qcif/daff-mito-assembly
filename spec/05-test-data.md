## 5. Test data

**Two-tier structure:**

| Tier | Trigger | Fixtures | What it validates |
|---|---|---|---|
| Fast CI | every push / PR | none (stub-only) | DSL syntax, channel wiring, container pulls, params load |
| Integration | nightly + on-demand | Tier 2 (Azure blob) | real tools, real assembly, biology assertions |

Per assembly target (`animal_mt`, `plant_pt`, `plant_mt`), integration
tests need one clean fixture. Contamination scenarios are added to
integration only when [P5](06-phases.md) starts. Skim depth
(~2× nuclear-equivalent) throughout — see [brief.md §2](brief.md).

| Assembly target | Fixture source | Notes |
|---|---|---|
| `animal_mt` | *Acyrthosiphon pisum* (pea aphid) SRR8306868 pre-recruited | Real ONT MinION WGS; ~2000 recruited reads → assembly-viable coverage. |
| `plant_pt` | *Datura stramonium* SRR11315861 pre-recruited | ~3000 recruited reads. |
| `plant_mt` | *Datura stramonium* SRR11315861 pre-recruited | ~1500 recruited reads. |

Acceptance criteria per integration fixture:
- Clean: organelle assembled to within target-appropriate length bounds,
  ≥ N expected loci extractable in `barcodes.fasta`, full annotation
  produced by `ANNOTATE`.
- Contamination (P5): target dominant assembly selected; any
  low-coverage secondary contigs recorded in diagnostics; target
  assembly unaffected by trace non-target reads.

See [task 11](../tasks/completed/11_integration_tests.md) for the integration
harness (Azure blob fixture layout, fetch script, assertions).

### 5.1 The integration fixtures are correctness tests, not benchmarks

**No number derived from these fixtures describes how well the workflow
performs.** They exist to prove the pipeline runs end to end, that stages
are wired correctly, and that decisions come out the way the spec says
they should. They are not a performance measurement and must never be
quoted as one. Four reasons, each independently disqualifying:

1. **They are deliberately small.** The whole run must fit a nightly
   GitHub Actions budget — 60 minutes on 2 vCPU / ~14 GB
   ([§5b](#5b-ci)). The fixtures are subsampled to hit that, so
   recruited-read counts are in the hundreds to low thousands. Real
   submissions are not.
2. **They are old public SRA data of unrecorded vintage.** The pipeline
   is *designed* for raw Dorado SUP R10.4.1 reads — that assumption is
   baked into the Flye read-mode default
   ([§9 item 5](07-open-questions.md#9-fine-tuning-post-prototype-benchmarking))
   and the CHOPPER/FILTLONG thresholds. The fixtures predate that
   chemistry and basecaller. We are validating a modern-ONT pipeline on
   pre-modern-ONT reads, and their instrument, chemistry, basecaller and
   date are **not currently recorded anywhere** — that should be fixed
   (see below).
3. **`n` = 1 per target, and it is really `n` = 2 samples, not 3.**
   `plant_pt` and `plant_mt` are the *same accession* (SRR11315861)
   pointed at two targets. There is no biological replication, no
   lineage diversity, and nothing that would reveal a threshold overfit
   to one genome.
4. **One of them is under-sequenced on its own declared target.**
   `INT-PLANT-01-mt` soft-fails the coverage gate and does not reach
   assembly at all ([`tasks/todo.md`](../tasks/todo.md)), so the entire
   `plant_mt` branch downstream of C2 has no integration coverage.

**Expect real performance to be substantially better.** Higher depth and
higher-quality reads should improve every metric these fixtures report —
recruitment yield, assembly contiguity, ORF integrity, locus recovery.
Fixture numbers are a **floor**, not an estimate.

**What this means in practice:**

- Every threshold currently in `nextflow.config` calibrated on these
  fixtures — `bin_target_thresholds`, `recruit_thresholds`,
  `coverage_limits` — is **provisional**, and is flagged as such at its
  definition. They are placeholders that happen to pass, not values
  anyone has validated.
- The fixture bounds in `tests/integration/expected/*/` pin *current
  observed behaviour* so regressions surface. They are not targets, and
  a change that moves a fixture outside its bounds is a finding to
  investigate rather than a bound to retune.
- **Threshold tuning against fixture bounds is circular and has already
  cost us.** Task 28 §10.2 could not raise `recruit_thresholds` above
  zero, not because zero is right, but because every non-zero value
  pushed a fixture out of its band — a band set by the same
  under-powered data. The knob is inert until better data exists.
- [§9](07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)'s
  benchmarks require the **P3 datasets**, not these fixtures. Where a §9
  item's experiment cannot be run on P3 data, the honest outcome is to
  say so and leave the item open — see
  [task 26](../tasks/26_binning_marker_genes.md), whose go/no-go gate is
  blocked for exactly this reason.

**Two things owed here:**

- **Record fixture provenance.** Each `tests/integration/expected/*/`
  directory should carry the accession's instrument, flowcell
  chemistry, basecaller + model, and run date alongside the organism and
  SRA ID it already records. Without that we cannot reason about how far
  fixture behaviour should generalise, or notice when it stops being
  representative.
- **Source a proper benchmark set.** Modern ONT (R10.4.1, SUP-basecalled)
  at realistic submission depth, with more than one taxon per target and
  `plant_pt` / `plant_mt` drawn from *different* samples. That set is
  what §9 tunes against; the Tier 2 fixtures stay as they are, doing the
  job they are actually good at.

## 5a. Tests

Three testing surfaces:

- **Python (`bin/*.py`) — `pytest`.** One test file per module under
  `scripts/tests/`, imported directly (no Nextflow harness). Each
  [custom-logic component](02-stages.md#22-custom-logic-components) (C1–C7)
  gets its own test module. Fixtures inline in the test file (via
  `tmp_path`). `flake8` runs alongside (per user CLAUDE.md). **Runs
  in fast CI on every push.**
- **Nextflow wiring (`main.nf`, `modules/local/*.nf`) — `-stub-run`.**
  Each process's `stub:` block emits `touch`ed outputs; the workflow
  runs end-to-end in ~60 s. Validates DSL syntax, channel connections,
  container pulls, and params loading. **Runs in fast CI on every
  push** ([task 10](../tasks/completed/10_ci_stub_only.md)).
- **Real-tool + biology — `-profile integration`.** Real tools run
  against Tier 2 fixtures fetched from Azure blob; assertions on
  assembly length, barcode recovery, annotation quality. **Runs
  nightly** ([task 11](../tasks/completed/11_integration_tests.md)).

**Rule of thumb for coverage:**

- Custom-logic (C1–C7): unit-tested with `pytest` at 100% branch
  coverage. These are our code; they must be exhaustively tested.
- Off-the-shelf tool wrappers (`NANOPLOT_*`, `CHOPPER`, `FILTLONG`,
  `RECRUIT`, `METAFLYE`, `BANDAGE_NG`, `BLAST_VALIDATE`,
  `ANNOTATE`, `ORGANELLE_MAP`, `MINIPROT_EXTRACT`): channel wiring
  covered by `-stub-run`; command-line correctness + output-shape
  handling covered by nightly integration.
- Network fetches (`scripts/build_refs.sh` calls to NCBI FTP): the
  refdata build script is out-of-band — not exercised by pipeline
  CI. Test at the script level with mocked HTTP if regression
  proofing is needed.

**No nf-test.** The `-stub-run` + integration split gives adequate
coverage without adding a third test framework. Revisit if per-process
regressions become common that neither tier catches.

## 5b. CI

CI runs on **GitHub Actions**. Two concerns: run every test surface
defined in [§5a](#5a-tests), and publish the bespoke container images
listed in [§2.2](02-stages.md#22-custom-logic-components) so the workflow has
something to `docker pull`.

**Workflows** (one YAML file each, under `.github/workflows/`):

- **`tests.yml`** — runs on every push and PR to any branch. Fast tier
  ([task 10](../tasks/completed/10_ci_stub_only.md)).
  1. `nextflow lint main.nf` — syntax + convention.
  2. **Container-coverage check** — assert every process in
     `modules/local/` resolves to a `container` directive via
     `conf/containers.config`, and that no tag is `latest`
     ([§1a](01-pipeline-flow.md#1a-engineering-constraints)). A process without a
     container fails CI.
  3. **Stageable-file check** — grep-based lint asserts no process
     `script:` block contains a bare `${params.<file>}` interpolation
     that isn't wrapped in `file(...)` or on the scalar-param allow-list
     ([§1a](01-pipeline-flow.md#1a-engineering-constraints)).
  4. **Python** — `pytest scripts/tests/` on `bin/*.py` +
     `flake8`, run inside the shared `neoformit/daff-wf5-scripts:test`
     image (`scripts/pytest.sh`), not host Python — C3's tests need
     `mappy` pinned to the same minimap2 version as the BIN_TARGET
     container, so one shared test image is simpler than per-component
     host deps.
  5. **Nextflow** — `nextflow run . -profile stub -stub-run`.
     End-to-end DSL + channel + container validation in ~60 s.
     Container runtime enabled so image-pull failures surface here.

- **`integration.yml`** — nightly + `workflow_dispatch`. Slow tier
  ([task 11](../tasks/completed/11_integration_tests.md)).
  1. Fetch Tier 2 fixtures from Azure blob.
  2. Fetch refdata bundle (cached).
  3. `nextflow run . -profile integration`.
  4. Run `tests/integration/assertions.sh` — biology assertions.

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
