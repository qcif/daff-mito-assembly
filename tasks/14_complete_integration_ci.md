# Task 14 — Complete the CI integration test infrastructure

**Phase:** P1 — closes out the integration-test scaffolding started in
[task 11](completed/11_integration_tests.md) and
[task 12](completed/12_azure_integration_fixtures.md).
**Goal:** Publish the reference bundle to Azure, wire up
`scripts/fetch_refs.sh`, and get `integration.yml` to actually run the
pipeline end-to-end (or fail at a genuine pipeline stage, not the fetch
step).

**Prerequisites:**
- [task 3](completed/3_refdata.md) — reference bundle exists locally at
  `refs/v2026.07/` (865 MB, `recruit/`, `validate/`, `proteins/`).
- [task 12](completed/12_azure_integration_fixtures.md) — Azure
  containers created with public blob access; fixtures fetch works.
- [task 13](13_fix_fast_ci.md) — fast CI green (independent, but the
  `-with-report`/`-with-trace` fix in §5 below applies the same lesson).

**Exit criteria:**

- Azure blob `refdata-wf5/v2026.07/refs.tar.gz` (single tarball) exists
  and is publicly downloadable.
- `scripts/fetch_refs.sh <version>` downloads and unpacks the bundle to
  `refs/<version>/`, verifies SHA256 against a committed manifest, and
  succeeds on a fresh checkout with no credentials.
- `.github/workflows/integration.yml` no longer references
  `AZURE_REFDATA_SAS_TOKEN`.
- A manual `gh workflow run integration.yml` reaches the **"Run
  pipeline"** step and either (a) completes all real-tool stages that
  are implemented, or (b) fails inside a specific Nextflow process with
  a diagnostic that points at the next unimplemented stage. It must
  **not** fail in `Fetch refdata bundle` or earlier.
- `docs/azure/setup.md` and `docs/testing.md` updated to reflect the
  credential-free workflow (drop the `AZURE_REFDATA_SAS_TOKEN` GitHub
  secret entirely; keep the batch start-task SAS token, which is
  separate).
- Progressive-assertion table below is committed as part of
  `tests/integration/assertions.sh` header comments so future stage
  tasks know which block to uncomment.

**Not in scope:**

- Building or re-building the reference bundle — [task 3](completed/3_refdata.md)
  already covers that.
- Bumping the bundle version (`v2026.07 → v2026.08`) — separate procedure
  when refseq / OGDraw datasets rev.
- Any real-tool wiring for stages still on stubs (RECRUIT, METAFLYE,
  ANNOTATE, VALIDATE, REPORT) — those are their own tasks.
- Caching the refdata bundle in GitHub Actions via `actions/cache`
  (marked as `# TODO` in `integration.yml` — deferred until we see the
  wall-clock cost of downloading 865 MB every run).

---

## 1. Package and publish the reference bundle

The bundle is large (~865 MB uncompressed). Package as a single
tarball rather than uploading each file individually — one HTTP
transfer, one SHA to verify, and it round-trips through GitHub
Actions cache cleanly if we add caching later.

```bash
cd /home/cameron/dev/daff-biosecurity/wf5

# Sanity-check bundle contents (should mirror plan.md §4.4 tree)
tree -L 2 refs/v2026.07 | head -40
du -sh refs/v2026.07

# Tar without absolute paths; deterministic sort for reproducible SHA
tar --sort=name \
    --owner=0 --group=0 --numeric-owner \
    --mtime='@0' \
    -czf refs-v2026.07.tar.gz \
    -C refs v2026.07

ls -lh refs-v2026.07.tar.gz
sha256sum refs-v2026.07.tar.gz
```

Upload to public blob (container already public-access blob per
task 12):

```bash
STORAGE_KEY=$(az storage account keys list \
    --account-name daffstandard \
    --resource-group daff-biosecurity \
    --query "[0].value" -o tsv)

az storage blob upload \
    --account-name daffstandard \
    --account-key "$STORAGE_KEY" \
    --container-name refdata-wf5 \
    --name "v2026.07/refs.tar.gz" \
    --file refs-v2026.07.tar.gz \
    --overwrite

# Verify anonymous download works
curl -sSI https://daffstandard.blob.core.windows.net/refdata-wf5/v2026.07/refs.tar.gz \
    | head -5
```

Commit the tarball SHA256:

```bash
sha256sum refs-v2026.07.tar.gz \
    | sed 's|refs-v2026.07.tar.gz|refs.tar.gz|' \
    > refs/v2026.07.sha256
git add refs/v2026.07.sha256
git commit -m "Pin SHA256 for refdata bundle v2026.07"
```

Then delete the local tarball (bundle stays uncompressed under `refs/`
for local development):

```bash
rm refs-v2026.07.tar.gz
```

## 2. Create `scripts/fetch_refs.sh`

Mirror `tests/integration/fetch_fixtures.sh` — plain `curl`, no
credentials, SHA verification.

```bash
#!/usr/bin/env bash
# Downloads and unpacks the reference bundle from Azure blob (public access).
#
# Usage:
#   bash scripts/fetch_refs.sh v2026.07

set -euo pipefail

VERSION="${1:?Usage: fetch_refs.sh <version>  e.g. v2026.07}"
STORAGE_ACCOUNT="${AZURE_STORAGE_ACCOUNT:-daffstandard}"
REFDATA_CONTAINER="${REFDATA_CONTAINER:-refdata-wf5}"
BASE_URL="https://${STORAGE_ACCOUNT}.blob.core.windows.net/${REFDATA_CONTAINER}"
TARBALL="refs-${VERSION}.tar.gz"

mkdir -p refs

if [ -d "refs/${VERSION}" ]; then
    echo "refs/${VERSION} already present — skipping download."
    exit 0
fi

echo "Downloading ${VERSION} refdata bundle..."
curl -fsSL "${BASE_URL}/${VERSION}/refs.tar.gz" -o "${TARBALL}"

echo "Verifying SHA256..."
sha256sum -c <(sed "s|refs.tar.gz|${TARBALL}|" "refs/${VERSION}.sha256")

echo "Unpacking..."
tar -xzf "${TARBALL}" -C refs
rm "${TARBALL}"

echo "Refdata ready in refs/${VERSION}/"
```

`chmod +x scripts/fetch_refs.sh`, commit.

## 3. Wire into `integration.yml`

Strip the `AZURE_REFDATA_SAS_TOKEN` env var (no longer needed) and
drop `-with-report`/`-with-trace` (they require `ps` inside every
process container — same lesson as [task 13](13_fix_fast_ci.md)):

```yaml
- name: Fetch refdata bundle
  # TODO: cache with actions/cache keyed on refs/v*.sha256
  run: bash scripts/fetch_refs.sh v2026.07

- name: Run pipeline
  run: nextflow run . -profile integration --outdir tests/integration/output
```

Also add a failure-only diagnostics upload (matching the pattern
proposed in [task 13](13_fix_fast_ci.md)):

```yaml
- name: Upload Nextflow failure diagnostics
  if: failure()
  uses: actions/upload-artifact@v4
  with:
    name: integration-nextflow-diagnostics
    path: |
      .nextflow.log
      work/**/.command.{sh,out,err,log,exitcode}
    if-no-files-found: ignore
    retention-days: 7
```

## 4. Documentation updates

**`docs/azure/setup.md`** — §5 (SAS tokens):

- Drop §5.1 "Refdata token (CI secret)" entirely.
- Renumber §5.2/5.3 to §5.1/5.2 (or just relabel — the batch start-task
  and pool JSON tokens are unchanged).
- Update §6 GitHub secrets table: `AZURE_REFDATA_SAS_TOKEN` no longer
  needed. If the batch start-task refdata token is also unnecessary now
  the container is public, drop §5.2 too and simplify.

**`docs/testing.md`** — §Tier 2 local run:

```bash
# One-off: fetch refdata bundle (public access, no credentials needed)
bash scripts/fetch_refs.sh v2026.07
```

Also drop the `AZURE_REFDATA_SAS_TOKEN` row from the required-secrets
table.

**`docs/azure/README.md`** — refresh any lingering SAS-token references
in the refdata rendering steps (§ "First-time setup" 3–6). The
`<REFDATA_SAS_TOKEN>` in `setup.sh.template` is separate (batch node
start task) — evaluate whether that also becomes public-URL-only.

## 5. Verify

```bash
# Smoke-test locally
rm -rf refs/v2026.07
bash scripts/fetch_refs.sh v2026.07
ls refs/v2026.07/recruit refs/v2026.07/validate refs/v2026.07/proteins | head

# Trigger CI
git push
gh workflow run integration.yml
gh run watch --exit-status
```

Expected result: the run reaches "Run pipeline". The pipeline itself
may fail inside a specific stage (e.g. RECRUIT if its module still
stubs a real command), but the failure diagnostic must clearly
identify the stage — not a missing script or credential.

## 6. Progressive assertion plan

Add this table as a header comment in
`tests/integration/assertions.sh` so the block-per-stage uncomment
sequence is discoverable from the file itself:

| Stage lands (task) | Uncomment block in `assertions.sh` |
|---|---|
| RECRUIT (real minimap2) | Coverage gate status per sample |
| METAFLYE (real assembler) | Assembly-length bounds vs `expected/*/assembly_bounds.json` |
| ANNOTATE (real annotator) | Barcode locus presence vs `expected/*/expected_loci.txt` |
| REPORT (real Jinja render) | `run-report.html` size > 100 KB, `metadata.json` schema |
| VALIDATE (real BLAST) | Taxonomic-identification consistency checks |

Each block in `assertions.sh` should include a `# TODO(task-N):
uncomment when X lands` line pointing at its owning task, so `grep TODO
tests/integration/assertions.sh` gives the outstanding work list.

## 7. Deliverables checklist

- [ ] `refs-v2026.07.tar.gz` uploaded to
      `daffstandard/refdata-wf5/v2026.07/refs.tar.gz`.
- [ ] `refs/v2026.07.sha256` committed.
- [ ] `scripts/fetch_refs.sh` committed and executable.
- [ ] `integration.yml` no longer references `AZURE_REFDATA_SAS_TOKEN`;
      `-with-report`/`-with-trace` removed; failure diagnostics step added.
- [ ] `docs/azure/setup.md`, `docs/testing.md`, `docs/azure/README.md`
      updated.
- [ ] Manual `integration.yml` run reaches "Run pipeline" and the
      failure (if any) is inside a pipeline process, not a fetch step.
- [ ] `AZURE_REFDATA_SAS_TOKEN` GitHub secret deleted
      (`gh secret delete AZURE_REFDATA_SAS_TOKEN`) once the workflow no
      longer references it.
- [ ] `assertions.sh` header updated with the progressive uncomment
      table from §6.
- [ ] Task moved to `tasks/completed/`.
