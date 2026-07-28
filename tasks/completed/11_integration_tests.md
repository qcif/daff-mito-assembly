# Task 11 — Integration test profile

**Phase:** P1 / P2 — the biology-validation counterpart to the
`-stub-run` fast-CI setup from [task 10](10_ci_stub_only.md).
**Goal:** Add `-profile integration` that runs the pipeline against
real Tier 2 fixtures (pre-recruited, assembly-viable) and validates
the full workflow end-to-end. Scheduled nightly; on-demand via
`workflow_dispatch`.

**Prerequisite:** [task 10](10_ci_stub_only.md) must land first
(deletes the tiny fixtures that would confuse the integration setup).

**Exit criteria:**

- `nextflow run main.nf -profile integration` completes end-to-end on
  a fresh checkout **after** `tests/integration/fetch_fixtures.sh`
  populates the fetched-fixture directory. Runtime target: ≤ 60 min
  wall-clock on GitHub Actions `ubuntu-latest` (4 vCPU, 16 GB).
- All 3 samples (`animal_mt`, `plant_pt`, `plant_mt`) produce
  non-empty organelle assemblies (see §4 assertions).
- Assertions defined in `tests/integration/assertions.sh` pass:
  assembly length in target-appropriate range, ≥ 3 expected
  loci in barcodes.fasta, valid metadata.json.
- Nightly GitHub Actions workflow (`.github/workflows/integration.yml`)
  green.
- No churn to fast CI (`.github/workflows/tests.yml` from
  [task 10](10_ci_stub_only.md)) — integration is additive.

**Not in scope:**

- Downstream-stage implementations (METAFLYE, MEDAKA, BIN_TARGET,
  BLAST_VALIDATE, ANNOTATE, MINIPROT_EXTRACT, ORGANELLE_MAP, COLLATE,
  RUN_REPORT). Each lands in its own task; this task sets up the
  harness so those tasks have somewhere to run assertions. Initial
  integration assertions therefore only cover what's real at the time
  task 11 lands (P1: stages 1-6 through COVERAGE_GATE).
- The production refdata bundle (published by [task 3](3_refdata.md)).
  Integration tests use the same `params.organelle_refs` /
  `params.blast_db` / `params.locus_panel` bundle — this task doesn't
  rebuild it.
- Azure infrastructure setup — reuses existing `deploy/azure/`
  machinery + SAS token flow from
  [`docs/azure/README.md`](../docs/azure/README.md).

**Cross-cutting rules:**

- Fixture reproducibility: every fetched FASTQ has a pinned SHA256
  in the fetch script; a mismatch is a fetch failure, not a warning.
- Assertions are shell-scriptable and independent of the Nextflow
  harness — same "local re-run" rule as fast CI ([spec §5b](../spec/05-test-data.md)).

---

## 1. Fixture layout

```
tests/integration/
├── samples.csv                     # committed — 3 rows, real read filenames
├── fetch_fixtures.sh               # committed — downloads fetched/ from Azure
├── expected/                       # committed — per-sample expected shapes
│   ├── animal_mt/
│   │   ├── expected_loci.txt       # COX1, CYTB, ND1, …
│   │   └── assembly_bounds.json    # { min_bp: 14000, max_bp: 20000 }
│   ├── plant_pt/
│   │   └── …
│   └── plant_mt/
│       └── …
├── assertions.sh                   # committed — post-run validation
├── staging/                        # gitignored — raw SRA (moved from tests/data/staging/)
│   └── README.md                   # committed — SRA metadata (was tests/data/staging/README.md)
└── fetched/                        # gitignored — Tier 2 fixtures from Azure
    ├── animal_mt.fastq.gz          # ~10-20 MB, pre-recruited
    ├── plant_pt.fastq.gz           # ~20-40 MB
    └── plant_mt.fastq.gz           # ~10-30 MB
```

`samples.csv`:

```csv
sample_id,assembly_target,reads,sample_info,sample_type,sample_receipt_date,storage_location
INT-ANIMAL-01,animal_mt,animal_mt.fastq.gz,Acyrthosiphon pisum SRR8306868 pre-recruited,intact,2019-05-02,ENA
INT-PLANT-01-pt,plant_pt,plant_pt.fastq.gz,Datura stramonium SRR11315861 pre-recruited,leaf,2021-04-02,ENA
INT-PLANT-01-mt,plant_mt,plant_mt.fastq.gz,Datura stramonium SRR11315861 pre-recruited,leaf,2021-04-02,ENA
```

Add to `.gitignore`:

```
tests/integration/fetched/
tests/integration/staging/*
!tests/integration/staging/README.md
```

## 2. `conf/integration.config`

```groovy
// Integration profile — real Tier 2 fixtures + production refs.
// Runs nightly (see .github/workflows/integration.yml).

params {
    samplesheet    = "${projectDir}/tests/integration/samples.csv"
    data_dir       = "${projectDir}/tests/integration/fetched"
    organelle_refs = "${projectDir}/refs/v2026.07/recruit"
    locus_panel    = "${projectDir}/refs/v2026.07/proteins"
    blast_db       = "${projectDir}/refs/v2026.07/validate"
    outdir         = "${projectDir}/tests/integration/output"
}

docker {
    enabled    = true
    runOptions = '-u $(id -u):$(id -g) --platform linux/amd64'
}
```

Refs paths assume [task 3](3_refdata.md) has published the bundle to
`refs/v2026.07/`. The integration workflow's fetch step must also
fetch the refs bundle (or use the same `fetch_fixtures.sh` to do
both). Consider: refs are large (~15 GB); fetching per run is
wasteful; cache in Azure blob or use `actions/cache`.

Register the profile in `nextflow.config`:

```groovy
profiles {
    stub                 { includeConfig 'conf/stub.config' }
    integration          { includeConfig 'conf/integration.config' }
    azure                { includeConfig 'conf/azure.config' }
    docker               { … }
    singularity          { … }
}
```

## 3. `tests/integration/fetch_fixtures.sh`

Downloads fetched fixtures + optionally the refs bundle from Azure blob
storage. Uses `azcopy` inside a container to avoid host tool dependency.

```bash
#!/usr/bin/env bash
# Downloads Tier 2 integration fixtures from Azure blob storage.
# Requires env var AZURE_FIXTURES_SAS_URL (set in .env.azure or CI secret).

set -euo pipefail

FIXTURES_DIR="tests/integration/fetched"
mkdir -p "$FIXTURES_DIR"

# One-shot azcopy inside a container
docker run --rm -u $(id -u):$(id -g) \
    -v "$PWD/$FIXTURES_DIR":/out \
    mcr.microsoft.com/azure-cli:2.60.0 \
    az storage blob download-batch \
        --sas-token "${AZURE_FIXTURES_SAS_TOKEN}" \
        --source integration-fixtures/v2026.07/ \
        --destination /out

# Verify SHA256s against pinned manifest
python3 tests/integration/verify_shas.py \
    "$FIXTURES_DIR" \
    tests/integration/fetched.sha256
```

`tests/integration/fetched.sha256` — committed manifest (small):

```
<sha256>  animal_mt.fastq.gz
<sha256>  plant_pt.fastq.gz
<sha256>  plant_mt.fastq.gz
```

Regenerate this file whenever fixtures change; SHA mismatch fails the
integration run.

## 4. `tests/integration/assertions.sh`

Runs after the Nextflow invocation. Fails the workflow on any
assertion miss. Assertions are progressive — start with existence
checks in P1, add biology checks as downstream stages ship.

```bash
#!/usr/bin/env bash
set -euo pipefail

OUTDIR="tests/integration/output"
FAILED=0

# --- Per-sample structural checks (always applicable) ------------
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    for f in metadata.json report.html; do
        if [[ ! -s "$OUTDIR/$sample/$f" ]]; then
            echo "FAIL: $sample/$f missing or empty"
            FAILED=1
        fi
    done
done

# --- Run-level checks -------------------------------------------
for f in run_manifest.json run-report.html; do
    if [[ ! -s "$OUTDIR/$f" ]]; then
        echo "FAIL: $f missing or empty"
        FAILED=1
    fi
done

# --- Biology checks (guarded by stage-implemented flag) ----------
# Uncomment as METAFLYE, MINIPROT_EXTRACT, ANNOTATE land.
#
# for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
#     asm="$OUTDIR/$sample/organelle_assembly.fasta"
#     if [[ ! -s "$asm" ]]; then
#         echo "FAIL: $sample assembly missing"
#         FAILED=1
#         continue
#     fi
#
#     asm_bp=$(grep -v '^>' "$asm" | tr -d '\n' | wc -c)
#     bounds="tests/integration/expected/${sample#INT-*-}/assembly_bounds.json"
#     min_bp=$(jq .min_bp "$bounds")
#     max_bp=$(jq .max_bp "$bounds")
#     if (( asm_bp < min_bp || asm_bp > max_bp )); then
#         echo "FAIL: $sample assembly $asm_bp bp outside [$min_bp, $max_bp]"
#         FAILED=1
#     fi
#
#     barcodes="$OUTDIR/$sample/barcodes.fasta"
#     expected="tests/integration/expected/${sample#INT-*-}/expected_loci.txt"
#     for locus in $(cat "$expected"); do
#         if ! grep -q ">${locus}" "$barcodes"; then
#             echo "WARN: $sample missing barcode $locus"
#             # Not a hard fail — locus recovery is best-effort
#         fi
#     done
# done

exit $FAILED
```

## 5. `.github/workflows/integration.yml`

Runs nightly at 02:00 UTC + on-demand via `workflow_dispatch` (button
in GitHub UI) + on any push to `main` that touches
`main.nf`/`modules/local/`/`conf/`.

```yaml
name: Integration tests

on:
  schedule:
    - cron: '0 2 * * *'
  workflow_dispatch:
  push:
    branches: [main]
    paths:
      - 'main.nf'
      - 'modules/local/**'
      - 'conf/**'
      - 'bin/**'
      - 'tests/integration/**'

jobs:
  integration:
    runs-on: ubuntu-latest
    timeout-minutes: 90
    steps:
      - uses: actions/checkout@v4

      - name: Install Nextflow
        uses: nf-core/setup-nextflow@v2
        with:
          version: '25.04.0'

      - name: Fetch integration fixtures
        env:
          AZURE_FIXTURES_SAS_TOKEN: ${{ secrets.AZURE_FIXTURES_SAS_TOKEN }}
        run: bash tests/integration/fetch_fixtures.sh

      - name: Fetch refdata bundle
        # TODO: cache with actions/cache keyed on refs bundle version
        env:
          AZURE_REFDATA_SAS_TOKEN: ${{ secrets.AZURE_REFDATA_SAS_TOKEN }}
        run: bash scripts/fetch_refs.sh v2026.07

      - name: Run pipeline
        run: |
          nextflow run . -profile integration \
              -with-report tests/integration/output/nextflow_report.html \
              -with-trace tests/integration/output/nextflow_trace.txt

      - name: Assertions
        run: bash tests/integration/assertions.sh

      - name: Upload outputs
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: integration-outputs
          path: tests/integration/output/
          retention-days: 7
```

## 6. Fixture generation recipe

Tier 2 fixtures are published to Azure blob storage under
`integration-fixtures/v<YYYY.MM>/` once — not regenerated per run.
Recipe below documents how to generate them; run manually when the
fixtures need refreshing.

```bash
STAGING="tests/integration/staging"  # raw SRA files, already downloaded
REFS="refs/v2026.07/recruit"          # production refs bundle
OUT="tests/integration/fetched"

RECRUIT_IMG=quay.io/biocontainers/mulled-v2-…
SEQTK_IMG=quay.io/biocontainers/seqtk:1.4--he4a0461_2

# Recruit + subsample per target. Target size: 10-40 MB per fixture.
# Read counts tuned to hit that size (approximately 500-3000 reads).

for target in animal_mt plant_pt plant_mt; do
    case $target in
        animal_mt) SRC=$STAGING/SRR8306868_1.fastq.gz; N=2000 ;;
        plant_pt)  SRC=$STAGING/SRR11315861_1.fastq.gz; N=3000 ;;
        plant_mt)  SRC=$STAGING/SRR11315861_1.fastq.gz; N=1500 ;;
    esac

    docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
        -v "$SRC":/reads.fastq.gz:ro \
        -v "$PWD/$REFS":/refs:ro \
        -v "$PWD/$OUT":/out \
        $RECRUIT_IMG bash -c "
            minimap2 -ax map-ont -t 4 /refs/${target}.mmi /reads.fastq.gz \\
                | samtools view -F 4 -q 1 -@ 4 \\
                | cut -f1 | sort -u > /out/${target}_ids.txt
            seqtk subseq /reads.fastq.gz /out/${target}_ids.txt \\
                | gzip > /out/${target}_recruited.fastq.gz"

    docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
        -v "$PWD/$OUT/${target}_recruited.fastq.gz":/in.fastq.gz:ro \
        -v "$PWD/$OUT":/out \
        $SEQTK_IMG bash -c "
            seqtk sample -s42 /in.fastq.gz $N | gzip > /out/${target}.fastq.gz"

    rm "$OUT/${target}_recruited.fastq.gz" "$OUT/${target}_ids.txt"
done

# Refresh sha256s
sha256sum $OUT/*.fastq.gz > tests/integration/fetched.sha256

# Publish to Azure blob
az storage blob upload-batch \
    --sas-token "$AZURE_FIXTURES_PUBLISH_SAS_TOKEN" \
    --destination integration-fixtures/v2026.07/ \
    --source $OUT
```

Size targets tuned to hit assembly-viable coverage:
- `animal_mt`: ~2000 recruited reads × 5 kb ≈ 10 MB total, ~600× mt coverage
- `plant_pt`: ~3000 × 15 kb ≈ 45 MB, ~300× plastid coverage
- `plant_mt`: ~1500 × 15 kb ≈ 22 MB, ~55× plant-mt coverage

Adjust `N` per target if fixture size targets shift.

## 7. Migration from `tests/data/staging/`

Move (with `git mv` for tracked files, plain `mv` for git-ignored):

```bash
git mv tests/data/staging/README.md tests/integration/staging/README.md
mv tests/data/staging/*.fastq.gz tests/integration/staging/
rmdir tests/data/staging tests/data  # empty after task 10 cleanup
```

Update `.gitignore` — the pattern was:
```
tests/data/staging/*
!tests/data/staging/README.md
```

Change to:
```
tests/integration/staging/*
!tests/integration/staging/README.md
tests/integration/fetched/
tests/integration/output/
```

## 8. Progressive assertion strategy

The `assertions.sh` biology block is commented-out on land. Uncomment
as stages become real:

| Stage completes | Uncomment |
|---|---|
| [task 8](8_recruit.md) RECRUIT (done) | Recruited FASTQ per sample non-empty |
| task 10-COV COVERAGE_GATE | `sample_status.json` `status == "ok"` on integration fixtures |
| task METAFLYE | `organelle_assembly.fasta` exists, length in expected range |
| task ANNOTATE | `organelle_annotation.gff` non-empty, contains expected gene features |
| task MINIPROT_EXTRACT | `barcodes.fasta` contains ≥ N expected loci |
| task COLLATE | `metadata.json` schema-valid, contains all upstream stage outputs |
| task RUN_REPORT | `run-report.html` renders + summary table has all 3 samples |

Each per-stage task's "Deliverables checklist" should include: "add
integration assertion for this stage in `tests/integration/assertions.sh`".

## 9. Deliverables checklist

- [ ] `conf/integration.config` — new profile.
- [ ] `tests/integration/samples.csv` — 3 rows (animal_mt, plant_pt, plant_mt).
- [ ] `tests/integration/fetch_fixtures.sh` — Azure blob download.
- [ ] `tests/integration/fetched.sha256` — pinned digests.
- [ ] `tests/integration/expected/` — per-sample expected shapes
      (assembly bounds, loci lists).
- [ ] `tests/integration/assertions.sh` — structural checks (P1);
      biology checks stubbed / commented for progressive uncomment.
- [ ] `tests/integration/staging/README.md` — moved from `tests/data/staging/`.
- [ ] `.gitignore` — updated for new paths.
- [ ] `nextflow.config` — `integration` profile registered.
- [ ] `.github/workflows/integration.yml` — nightly + workflow_dispatch.
- [ ] Azure blob container `integration-fixtures/v2026.07/` populated
      (one-off, see §6 recipe).
- [ ] GitHub secret `AZURE_FIXTURES_SAS_TOKEN` set.
- [ ] `docs/testing.md` (new) — two-tier CI overview.
- [ ] First scheduled run green.
- [ ] `tasks/completed.txt` appended.

## 10. Notes / non-issues

- **Refs bundle fetch is the slow step, not the fixture fetch.** ~15 GB
  refs vs ~80 MB fixtures. Cache aggressively; consider mounting from
  Azure blob directly (nextflow-azure has FUSE support) rather than
  downloading. Iterate on this after the first run to see where the
  bottleneck lands in practice.
- **Assertions must be independent of the report renderer.** The
  report is a downstream consumer; assertions consume `metadata.json`
  directly. If the report fails to render, we still want per-stage
  assertions to fail-loudly on the actual defect.
- **Nightly-only cadence is a trade-off.** A developer landing a
  METAFLYE fix has to wait up to 24 h to know if it works on real
  data. `workflow_dispatch` covers on-demand cases; if it becomes a
  friction point, tighten the schedule (every 6 h) or add a
  `push:paths:['modules/local/metaflye.nf', …]` trigger for the
  hot-path files.
- **No `-profile integration_bad_samplesheet`.** Preflight failure
  coverage is exhaustive at the unit-test level
  ([`test_parse_samplesheet.py`](../scripts/tests/test_parse_samplesheet.py)),
  no value in re-testing at integration scale.
