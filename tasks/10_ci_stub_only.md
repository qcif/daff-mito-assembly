# Task 10 — Drop `-profile test`; CI becomes `-stub-run` only

**Phase:** P1 CI simplification.
**Goal:** Reduce test-fixture maintenance burden. Delete the tiny
pre-recruited real-tool fixtures and the `-profile test` /
`-profile test_bad_samplesheet` profiles that used them. Fast CI on
every push becomes `-stub-run` only (validates channel wiring +
container pulls + params); real-tool + biology validation moves to
[task 11 — integration tests](11_integration_tests.md), run nightly.

**Motivation** — captured in [task 8 §9](8_recruit.md#9-findings--recruitment-yield-on-ci-fixtures):
- With 70-read WGS fixtures, RECRUIT yields 18-86 organelle reads per
  sample. That's below any realistic COVERAGE_GATE MIN threshold, so
  once [task 10-COV — COVERAGE_GATE](.) lands (unwritten), every CI
  sample soft-fails and skips METAFLYE → BIN_TARGET → ORGANELLE_MAP
  via the failed-branch shortcut in [main.nf:104-108](../main.nf#L104-L108).
- METAFLYE on 56 reads produces nothing assemblable regardless of gate
  threshold, so the assembly-side wiring can never be exercised at CI
  fixture sizes.
- The fixture-generation pipeline (recruit + subsample from staging
  SRA files) is elaborate — regenerating it every time a threshold
  changes is real toil.
- Fast CI value with tiny fixtures diminishes to "chopper 0.14 changed
  its flag names" class of regressions. That is real value but not
  worth an ongoing fixture-maintenance burden when integration tests
  will catch the same thing within 24h.

**Exit criteria:**

- `nextflow run main.nf -profile stub` (new profile — see §1) completes
  end-to-end in ≤ 60 seconds. All processes hit their `stub:` blocks;
  every expected output file (per [spec §0](../spec/00-overview.md))
  is `touch`ed under `results/STUB-01/` (the single sample emitted by
  VALIDATE_SAMPLESHEET's stub).
- `-profile test` and `-profile test_bad_samplesheet` are removed;
  `nextflow run main.nf -profile test` fails with "unknown profile".
- All `tests/data/*.fastq.gz` and `tests/data/refs/` contents (except
  `tests/data/staging/README.md`) are deleted from git.
- `.github/workflows/tests.yml` is updated: the Nextflow job runs
  `-profile stub`, and the output-layout check asserts `STUB-01`
  (not `TEST-*`).
- Python unit tests ([`scripts/tests/test_parse_samplesheet.py`](../scripts/tests/test_parse_samplesheet.py))
  continue to pass — they run outside Nextflow and are unaffected.
- Total git delta from this task is a **reduction**: ~5-10 MB removed
  (fixtures + subset refs), ~50 lines of config removed.

**Not in scope:**

- Building the integration profile (nightly real-tool tests) — that is
  [task 11](11_integration_tests.md).
- Changing what `-stub-run` validates within Nextflow; every process
  already has a `stub:` block from [task 1 scaffold](1_scaffold.md).
- Removing `tests/data/staging/` — its `README.md` records the SRA
  accessions we'll reuse for integration fixtures. Staging FASTQs are
  already `.gitignore`'d; leave them on disk for whoever regenerates
  Tier 2 fixtures next.

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pulls still happen under `-stub-run` — stubs run inside
  their declared containers, so container-tag issues surface here.
- Unit tests remain the primary regression surface for
  [C1–C7 custom logic](../spec/02-stages.md#22-custom-logic-components).

---

## 1. New `-profile stub`

Replaces `-profile test`. Minimal fixture footprint — just enough to
satisfy Nextflow's `path`-input staging requirements. `-stub-run` is
implied by the profile so callers don't need to remember to pass it.

New file: `conf/stub.config`

```groovy
// Stub-only profile for fast CI wiring validation.
// Every process runs its stub: block; no real tools are invoked.
// Real-tool + biology validation lives in -profile integration
// (see tasks/11_integration_tests.md).

params {
    samplesheet    = "${projectDir}/tests/stub/samplesheet.csv"
    data_dir       = "${projectDir}/tests/stub/data"
    organelle_refs = "${projectDir}/tests/stub/refs/organelle_refs"
    locus_panel    = "${projectDir}/tests/stub/refs/locus_panel.faa"
    blast_db       = "${projectDir}/tests/stub/refs/blast_db"
    outdir         = "${projectDir}/tests/output"
}

process {
    resourceLimits = [ cpus: 2, memory: '4.GB', time: '10.m' ]
}

docker {
    enabled    = true
    runOptions = '-u $(id -u):$(id -g) --platform linux/amd64'
}
```

Callers invoke: `nextflow run main.nf -profile stub -stub-run`. The
`-stub-run` flag is deliberately explicit rather than baked into the
profile — Nextflow doesn't support `stub-run = true` as a config
directive, and hiding it behind a profile alias would obscure what's
happening.

### 1.1 New fixture layout: `tests/stub/`

Minimal placeholders that satisfy `path` input staging. Content is
never read.

```
tests/stub/
├── samplesheet.csv                 # 20 bytes — header row only
├── data/
│   └── .gitkeep
└── refs/
    ├── organelle_refs/
    │   └── .gitkeep
    ├── blast_db/
    │   └── .gitkeep
    └── locus_panel.faa             # 0 bytes
```

`samplesheet.csv` contains just the header row — no data rows needed
because VALIDATE_SAMPLESHEET's stub emits a hardcoded `STUB-01` sample
regardless of input:

```csv
sample_id,assembly_target,reads,sample_info,sample_type,sample_receipt_date,storage_location
```

## 2. Files to delete

Delete (with `git rm`):

- `conf/test.config`
- `conf/test_bad_samplesheet.config`
- `tests/data/samples.csv`
- `tests/data/samples_bad.csv`
- `tests/data/README.md`
- `tests/data/animal.fastq.gz`
- `tests/data/animal_b.fastq.gz`
- `tests/data/plant.fastq.gz`
- `tests/data/refs/recruit/` (whole directory — indices + FASTAs + README)
- `tests/data/refs/organelle_refs.mmi` (empty placeholder)
- `tests/data/refs/blast_db/` (empty placeholder — moved to `tests/stub/`)
- `tests/data/refs/locus_panel.faa` (empty placeholder — moved to `tests/stub/`)

**Keep**:
- `tests/data/staging/README.md` — records SRA accessions for
  reuse by [task 11](11_integration_tests.md). Move to
  `tests/integration/staging/README.md` once task 11 lands.
- `tests/data/staging/*.fastq.gz` — already `.gitignore`'d; leave on
  disk for local integration-fixture regeneration.

After cleanup the `tests/data/` tree collapses to just `staging/` and
that becomes the migration target for [task 11's](11_integration_tests.md)
integration harness.

## 3. `nextflow.config` — profile changes

Replace the `test` and `test_bad_samplesheet` profile registrations
with a single `stub` profile:

```groovy
profiles {
    stub                 { includeConfig 'conf/stub.config' }
    azure                { includeConfig 'conf/azure.config' }
    docker               { … }        // unchanged
    singularity          { … }        // unchanged
    // -profile integration lands via task 11
}
```

Do not add a `test` profile alias for `stub` — a hard error on the
old name surfaces stale docs / CI configs immediately.

## 4. `.github/workflows/tests.yml` updates

Existing workflow uses `-profile test -stub` (see
[current tests.yml](../.github/workflows/tests.yml)). Two changes:

1. **Nextflow-test job invocation** — swap:
   ```yaml
   - name: Run stub pipeline
     run: |
       ~/.local/bin/nextflow run . -profile stub -stub-run \
           --outdir tests/output \
           -with-report tests/output/nextflow_report.html \
           -with-trace tests/output/nextflow_trace.txt
   ```

2. **Output-layout check** — replace hardcoded `TEST-PLANT-01` /
   `TEST-ANIMAL-01` sample IDs with the single `STUB-01`:
   ```python
   for sample_id in ['STUB-01']:
       sample_dir = outdir / sample_id
       for fname in ['metadata.json', 'report.html']:
           …
   ```

3. **Container pull** — the current `docker pull` line pulls a
   pre-P1 coreutils image. Update to pull the actual containers used
   by processes-in-stub-mode (which is now most of the biocontainers
   for QC/RECRUIT plus `python:3.12-slim` for the C1-C7 modules — see
   [`conf/containers.config`](../conf/containers.config)). Or drop
   the explicit pre-pull entirely and let Nextflow pull on demand;
   with `-stub-run` the pulls take ~1 min total.

4. **Python-tests job** — no changes needed; unit tests are unaffected.
   Consider adding a pytest invocation for
   [`scripts/tests/test_parse_samplesheet.py`](../scripts/tests/test_parse_samplesheet.py)
   if that isn't already covered by the existing `unittest discover`.

## 5. `spec/05-test-data.md` update

The current spec describes tests as if `-profile test` runs real
tools. Rewrite §5a to reflect the two-tier structure:

- **Fast CI (this task):** `-stub-run` end-to-end + pytest unit tests
  on `bin/*.py`. No real tools executed by the pipeline.
- **Nightly CI ([task 11](11_integration_tests.md)):** `-profile
  integration` with fetched Tier 2 fixtures. Real tools, full
  assembly, assertions on biology.

Update the "Rule of thumb for what to mock" section: the tool-wrapper
processes are no longer mock-tested at the Nextflow level in fast CI
— they run for real in integration only. Nextflow-level tests in
fast CI reduce to `-stub-run` end-to-end coverage.

## 6. Migration steps (execution order)

1. Create `conf/stub.config` and `tests/stub/` layout.
2. Update `nextflow.config` profiles block.
3. Update `.github/workflows/tests.yml`.
4. Update `spec/05-test-data.md`.
5. Delete files listed in §2.
6. Run `nextflow run main.nf -profile stub -stub-run` locally to
   confirm 62/62 → 17/17 processes complete (`STUB-01` × 17 stages).
7. Run `pytest scripts/tests/` — 20 passing (no change).
8. Commit + push; CI should be green.

## 7. Deliverables checklist

- [ ] `conf/stub.config` — new, minimal params.
- [ ] `tests/stub/{samplesheet.csv,data/.gitkeep,refs/…}` — new
      minimal fixtures.
- [ ] `nextflow.config` — profiles block updated.
- [ ] `.github/workflows/tests.yml` — invocation + output-layout
      assertions updated for `STUB-01`.
- [ ] `spec/05-test-data.md` — two-tier CI reflected.
- [ ] All files listed in §2 deleted via `git rm`.
- [ ] `tests/data/README.md` deleted (redundant once fixtures are
      gone); `tests/data/staging/README.md` moved into
      `tests/integration/staging/` by [task 11](11_integration_tests.md)
      or left in place until then.
- [ ] Local run: `nextflow run main.nf -profile stub -stub-run`
      completes; `results/STUB-01/{metadata.json,report.html}` exist.
- [ ] Local run: `pytest scripts/tests/` — 20 pass.
- [ ] CI green on push.
- [ ] `tasks/completed.txt` appended.

## 8. Notes / non-issues

- **VALIDATE_SAMPLESHEET stub content is authoritative.** With no real
  samplesheet content read in stub mode, the sample ID / target / meta
  fields the stub emits (`STUB-01`, `animal_mt`, empty everything else)
  are what downstream processes see. If the stub gets out of sync with
  the schema (added a required meta field, changed field names), stub
  CI breaks in a diagnostic way — that's a feature.
- **Container pulls still happen.** `-stub-run` executes each process's
  `stub:` block **inside** its container. So container-tag issues,
  registry outages, and platform mismatches still surface in fast CI.
  This is the main non-trivial thing stub CI still validates beyond
  what `nextflow lint` catches.
- **We lose the "duplicate sample_id fails preflight" end-to-end
  check.** That coverage moves entirely into
  [`test_parse_samplesheet.py::test_duplicate_sample_id_fails`](../scripts/tests/test_parse_samplesheet.py)
  (unit test). The Nextflow-level check was expensive to maintain for
  what it validated.
- **What we lose that unit tests can't replace:** verifying that
  VALIDATE_SAMPLESHEET's non-zero exit actually aborts the whole
  Nextflow run before downstream processes are scheduled. That's a
  Nextflow-behaviour claim, not a Python claim. Consider a targeted
  nf-test if this becomes flaky in practice — but not needed for
  this task.
