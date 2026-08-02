# Task 12 — Publish integration fixtures to Azure and wire CI secrets

**Phase:** P1 — prerequisite for the first nightly integration run.
**Goal:** Subsample the recruited staging files into Tier 2 fixtures,
upload them to Azure blob storage, compute pinned SHA256s, and make this
container public (anonymous blob access) so `.github/workflows/integration.yml`
can run end-to-end.

**Prerequisite:** [task 11](11_integration_tests.md) must
land first (`conf/integration.config`, `fetch_fixtures.sh`,
`fetched.sha256` placeholder all committed).

**Exit criteria:**

- `tests/integration/fetched/` populated locally with 3 fixture files
  (sizes within targets in §2).
- `tests/integration/fetched.sha256` updated with real SHA256s and
  committed.
- Azure blob container `test-data/wf5/v2026.07/` contains all
  3 fixture files (verified with `az storage blob list`).
- `bash tests/integration/fetch_fixtures.sh` succeeds on a fresh
  checkout (fixtures absent, blob present).
- `.github/workflows/integration.yml` first scheduled or
  `workflow_dispatch` run reaches the "Run pipeline" step (pipeline
  itself may still fail on missing refdata — that's expected until
  [task 3](3_refdata.md) publishes the production bundle).

**Not in scope:**

- Publishing the production refdata bundle to `refdata-wf5` — that
  is [task 3](3_refdata.md).
- Generating integration fixtures for contamination scenarios (P5).

---

## 1. Azure storage layout

Both containers live on the **`daffstandard`** storage account
(standard LRS, already used for `scripts/` and `workdata/`). Both should have
public access enabled so we don't have to worry about auth - all the data here
is public anyway.

| Container | Purpose |
|---|---|
| `test-data` | Tier 2 recruited FASTQ fixtures, versioned by subfolder |
| `refdata-wf5` | Production reference data bundle (populated by task 3) |

Fixture path inside `test-data`:

```
test-data/
└── wf5/
    └── v2026.07/
        ├── animal_mt.fastq.gz   (~7 MB)
        ├── plant_pt.fastq.gz    (~40 MB)
        └── plant_mt.fastq.gz    (~20 MB)
```

## 2. Subsampling recipe

The recruited files in `tests/integration/staging/` are the source.
`animal_mt` is used as-is (all 1242 reads); `plant_pt` and `plant_mt`
are subsampled with `seqtk sample` to hit size / coverage targets.

```bash
STAGING="tests/integration/staging"
OUT="tests/integration/fetched"
mkdir -p "$OUT"

SEQTK_IMG=quay.io/biocontainers/seqtk:1.4--he4a0461_2

# animal_mt — use all recruited reads (1242 reads, ~7 MB)
cp "$STAGING/animal_mt_recruited.fastq.gz" "$OUT/animal_mt.fastq.gz"

# plant_pt — subsample to 3000 reads (~40 MB)
docker run --rm -u "$(id -u):$(id -g)" --platform linux/amd64 \
    -v "$PWD/$STAGING/plant_pt_recruited.fastq.gz":/in.fastq.gz:ro \
    -v "$PWD/$OUT":/out \
    $SEQTK_IMG \
    bash -c 'seqtk sample -s42 /in.fastq.gz 3000 | gzip > /out/plant_pt.fastq.gz'

# plant_mt — subsample to 1500 reads (~20 MB)
docker run --rm -u "$(id -u):$(id -g)" --platform linux/amd64 \
    -v "$PWD/$STAGING/plant_mt_recruited.fastq.gz":/in.fastq.gz:ro \
    -v "$PWD/$OUT":/out \
    $SEQTK_IMG \
    bash -c 'seqtk sample -s42 /in.fastq.gz 1500 | gzip > /out/plant_mt.fastq.gz'
```

Verify sizes and read counts:

```bash
ls -lh tests/integration/fetched/
for f in animal_mt plant_pt plant_mt; do
    echo -n "$f: "
    zcat "tests/integration/fetched/${f}.fastq.gz" | awk 'NR%4==1' | wc -l
done
```

## 3. Compute and commit SHA256s

```bash
sha256sum tests/integration/fetched/*.fastq.gz \
    | sed 's|tests/integration/fetched/||' \
    > tests/integration/fetched.sha256

cat tests/integration/fetched.sha256
git add tests/integration/fetched.sha256
git commit -m "Add pinned SHA256s for Tier 2 integration fixtures v2026.07"
```

## 4. Upload fixtures to Azure blob

```bash
source deploy/azure/batch-helpers.sh   # loads helpers + env vars

# Create containers (idempotent) — public blob access so CI needs no credentials
az storage container create \
    --account-name daffstandard \
    --name test-data \
    --auth-mode login \
    --public-access blob

az storage container create \
    --account-name daffstandard \
    --name refdata-wf5 \
    --auth-mode login \
    --public-access blob

# Upload fixtures
for f in animal_mt plant_pt plant_mt; do
    az storage blob upload \
        --account-name daffstandard \
        --container-name test-data \
        --name "wf5/v2026.07/${f}.fastq.gz" \
        --file "tests/integration/fetched/${f}.fastq.gz" \
        --auth-mode login \
        --overwrite
done

# Verify
az storage blob list \
    --account-name daffstandard \
    --container-name test-data \
    --prefix wf5/v2026.07/ \
    --auth-mode login \
    --query "[].{name:name, size:properties.contentLength}" \
    -o table
```

## 7. Smoke-test the fetch script

On a fresh terminal (no local `tests/integration/fetched/` contents):

```bash
rm -rf tests/integration/fetched/
bash tests/integration/fetch_fixtures.sh
ls -lh tests/integration/fetched/
```

Expected: 3 files, SHA256 verification passes, script exits 0.

## 8. Trigger a manual integration run

```bash
gh workflow run integration.yml
gh run watch   # follow progress
```

The pipeline will proceed to "Run pipeline" and likely fail at the
`fetch_refs.sh` step (refdata bundle not yet published). That's
expected — the failure proves the fixture fetch and secret wiring work.
Assess the failure mode and open a follow-up if it fails earlier.

## 9. Deliverables checklist

- [ ] `tests/integration/fetched/` populated locally.
- [ ] `tests/integration/fetched.sha256` — real SHA256s, committed.
- [ ] Azure container `test-data` created.
- [ ] Azure container `refdata-wf5` created (empty; populated by task 3).
- [ ] Fixtures uploaded to `test-data/wf5/v2026.07/`.
- [ ] `fetch_fixtures.sh` smoke-test passes.
- [ ] Manual `integration.yml` dispatch reaches "Run pipeline" step.
- [ ] Task moved to `tasks/completed/`.
