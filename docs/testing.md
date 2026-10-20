# Testing

Two-tier CI: fast CI on every push, nightly integration on real data.

## Tier 1 — Fast CI (every push / PR)

**Trigger:** `.github/workflows/tests.yml`, runs on every push and PR.

**What it validates:**
- `nextflow lint main.nf` — DSL syntax and conventions.
- Container-coverage check — every process in `modules/local/` must resolve to a `container` directive in `conf/containers.config`. No `latest` tags.
- Bare-params file check — no `${params.<file>}` in `script:` blocks that isn't wrapped in `file()`.
- Python unit tests (`pytest scripts/tests/`) + `flake8` on `bin/` and `scripts/`.
- Nextflow stub end-to-end: `nextflow run . -profile stub -stub-run`. All 17 processes hit their `stub:` blocks, expected output files exist under `tests/output/STUB-01/`.

**Fixtures:** `tests/stub/` — minimal placeholders (header-only samplesheet, empty dirs). No real data read; VALIDATE_SAMPLESHEET's stub emits a hardcoded `STUB-01` sample.

**Runtime:** ~2 min.

**Run locally:**

```bash
nextflow run . -profile stub -stub-run
pytest scripts/tests/
flake8 scripts/ bin/ --max-line-length=100
```

## Tier 2 — Integration tests (nightly + on-demand)

**Trigger:** `.github/workflows/integration.yml`, runs at 02:00 UTC daily and on `workflow_dispatch`. Also runs on push to `main` touching `main.nf`, `modules/local/`, `conf/`, `bin/`, or `tests/integration/`.

**What it validates:**
- Real tools run against Tier 2 pre-recruited ONT fixtures.
- Full pipeline end-to-end (all 17 stages).
- Post-run assertions in `tests/integration/assertions.sh` — structural checks now; biology checks (assembly length, barcode locus recovery) uncommented progressively as downstream stages land.

**Fixtures:** `tests/integration/fetched/` — gitignored, fetched at CI run time from Azure blob storage. Three samples:

| Sample | Assembly target | Source | ~Size |
|---|---|---|---|
| `INT-ANIMAL-01` | `animal_mt` | *Acyrthosiphon pisum* SRR8306868 | ~10 MB |
| `INT-PLANT-01-pt` | `plant_pt` | *Datura stramonium* SRR11315861 | ~40 MB |
| `INT-PLANT-01-mt` | `plant_mt` | *Datura stramonium* SRR11315861 | ~20 MB |

**Runtime:** ≤ 60 min wall-clock (4 vCPU / 16 GB `ubuntu-latest`).

**Run locally** (requires Azure credentials and the refdata bundle):

```bash
# One-off: fetch Tier 2 fixtures from Azure blob (public access, no credentials needed)
bash tests/integration/fetch_fixtures.sh

# One-off: fetch refdata bundle (public access, no credentials needed)
bash scripts/fetch_refs.sh v2026.08

# Run integration profile
nextflow run . -profile integration

# Assert outputs
bash tests/integration/assertions.sh
```

## Fixture sources and regeneration

Tier 2 fixtures are pre-recruited subsamples of raw ONT WGS runs. Raw staging files live in `tests/integration/staging/` (gitignored). To regenerate:

1. Ensure staging files are present (see `tests/integration/staging/README.md` for ENA download URLs).
2. Recruit each target from its raw run the way the pipeline's RECRUIT stage does, against the current `refs/<version>/recruit` panel, then subsample. Run it in the `coverage-gate` image from `conf/containers.config` (minimap2, seqtk and python3):
   ```
   minimap2 -x map-ont -t 8 refs/<version>/recruit/<target>.mmi <raw>.fastq.gz > hits.paf
   python3 bin/recruit_filter.py --paf hits.paf --min-aligned-frac 0 \
       --min-aligned-bp 0 --stats recruit_stats.json --out ids.txt
   seqtk subseq <raw>.fastq.gz ids.txt | gzip > recruited.fastq.gz
   seqtk sample -s42 recruited.fastq.gz <N> | gzip > <target>.fastq.gz
   ```
   `<N>` is 2000 (`animal_mt`), 3000 (`plant_pt`) and 1500 (`plant_mt`); a target with fewer recruited reads than `<N>` keeps them all. Rebuild a fixture whenever the recruit panel for its target changes, so the fixture is not recruited against a different panel from the one the pipeline uses.
3. Upload all three to Azure blob under `test-data/wf5/v<YYYY.MM>/` (a new version directory holds the full set, so copy unchanged fixtures across) and set `FIXTURES_VERSION` in `tests/integration/fetch_fixtures.sh` and `DEFAULT_FIXTURES_VERSION` in `deploy/azure/run-wf5-azure.sh`.
4. Recompute SHA256s: `sha256sum tests/integration/fetched/*.fastq.gz > tests/integration/fetched.sha256`.
5. Commit the updated `fetched.sha256`.

## Progressive assertion strategy

Biology assertions in `tests/integration/assertions.sh` are commented out until the relevant stage is real. See `tasks/completed/11_integration_tests.md §8` for the per-stage uncomment plan.

## Required secrets (GitHub)

None. Both the integration fixtures and the refdata bundle are served from public Azure blob
containers and require no credentials.
