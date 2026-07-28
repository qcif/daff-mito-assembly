# Test fixtures

Real ONT reads for CI wiring tests. Small enough to commit; not large
enough to complete an actual organelle assembly — that's a Tier 2
integration concern, tracked separately.

## Files

All fixtures are **pre-recruited** — reads were first aligned against the
target organelle references (see `tests/data/refs/recruit/`), and the
fixtures are subsampled from those recruited reads. This guarantees
non-zero RECRUIT output in CI without needing large fixture files.

| File | SHA256 | Reads | Size | Purpose |
|---|---|---|---|---|
| `animal.fastq.gz` | `a588b63d63d368ad93cb4e75bce7bb49a0870332ddf012e430834c859bc9027b` | 70 | 380 K | `TEST-ANIMAL-01` primary; `TEST-ANIMAL-02` first file. Pre-recruited from `animal_mt` pool (1,242 reads available). |
| `animal_b.fastq.gz` | `dba2778d4fe308f2c3ba347095ed48cf2206f55f81f8ae064a71874cf2f62698` | 35 | 199 K | `TEST-ANIMAL-02` second file — seed 137 to avoid read-name collisions with `animal.fastq.gz`. |
| `plant.fastq.gz` | `673d61ab11afb7acd42890e46c83f483e290fe57c0fc89f3f24c181c196cc421` | 30 | 552 K | `TEST-PLANT-01-pt` and `TEST-PLANT-01-mt`. Mix of 15 reads from `plant_pt` pool (53,641 available) + 15 reads from `plant_mt` pool (14,981 available), ensuring non-zero recruitment for both targets. |

## Sources

- **Animal fixtures** ← `SRR8306868` — *Acyrthosiphon pisum* (pea
  aphid), ONT MinION WGS. 892,247 reads / 7.3 Gb raw. Published
  2019-05-02.
- **Plant fixture** ← `SRR11315861` — *Datura stramonium*
  (Jimsonweed), ONT MinION WGS. 558,821 reads / 8.4 Gb raw. Published
  2021-04-02.

Raw source files are not committed; they live under
[`tests/data/staging/`](staging/) which is git-ignored (see
[`.gitignore`](../../.gitignore)) except for
[`staging/README.md`](staging/README.md) which records metadata + ENA
download URLs.

## Regenerating the fixtures

Fixtures are pre-recruited: first recruit from the full SRA against the
test organelle refs, then subsample. Build the test refs first if
they don't exist (see `tests/data/refs/recruit/README.md`).

```bash
RECRUIT_IMG="quay.io/biocontainers/mulled-v2-9278ceb357570ba6e25522f5a16e2a4d3ba61a68:74b9019c6ae38f81f95dd09e981151bb2d1028ee-0"
SEQTK_IMG="quay.io/biocontainers/seqtk:1.4--he4a0461_2"
REFS="$PWD/tests/data/refs/recruit"
STAGING="$PWD/tests/data/staging"

# Step 1: Recruit from full SRA files into staging/
docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v "$STAGING/SRR8306868_1.fastq.gz":/reads.fastq.gz:ro \
    -v "$REFS":/refs:ro -v "$STAGING":/out \
    $RECRUIT_IMG bash -c '
        minimap2 -ax map-ont -t 4 /refs/animal_mt.mmi /reads.fastq.gz \
            | samtools view -F 4 -q 1 -@ 4 | cut -f1 | sort -u > /out/animal_mt_ids.txt
        seqtk subseq /reads.fastq.gz /out/animal_mt_ids.txt | gzip > /out/animal_mt_recruited.fastq.gz'

docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v "$STAGING/SRR11315861_1.fastq.gz":/reads.fastq.gz:ro \
    -v "$REFS":/refs:ro -v "$STAGING":/out \
    $RECRUIT_IMG bash -c '
        for target in plant_pt plant_mt; do
            minimap2 -ax map-ont -t 4 /refs/${target}.mmi /reads.fastq.gz \
                | samtools view -F 4 -q 1 -@ 4 | cut -f1 | sort -u > /out/${target}_ids.txt
            seqtk subseq /reads.fastq.gz /out/${target}_ids.txt | gzip > /out/${target}_recruited.fastq.gz
        done'

# Step 2: Subsample from recruited pools
docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v "$STAGING/animal_mt_recruited.fastq.gz":/in.fastq.gz:ro \
    -v "$PWD/tests/data":/out $SEQTK_IMG \
    bash -c 'seqtk sample -s42 /in.fastq.gz 70 | gzip > /out/animal.fastq.gz'

docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v "$STAGING/animal_mt_recruited.fastq.gz":/in.fastq.gz:ro \
    -v "$PWD/tests/data":/out $SEQTK_IMG \
    bash -c 'seqtk sample -s137 /in.fastq.gz 35 | gzip > /out/animal_b.fastq.gz'

docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v "$STAGING/plant_pt_recruited.fastq.gz":/pt.fastq.gz:ro \
    -v "$STAGING/plant_mt_recruited.fastq.gz":/mt.fastq.gz:ro \
    -v "$PWD/tests/data":/out $SEQTK_IMG \
    bash -c '
        seqtk sample -s42 /pt.fastq.gz 15 > /tmp/pt.fastq
        seqtk sample -s42 /mt.fastq.gz 15 > /tmp/mt.fastq
        cat /tmp/pt.fastq /tmp/mt.fastq | gzip > /out/plant.fastq.gz'
```

Seeds are pinned so anyone with the same source file gets the same
fixture. If read counts need to change (e.g. to hit a different size
target), pick a new count and update this file's SHA256 table.

## Sample sheet

[`samples.csv`](samples.csv) references these fixtures. Per
[spec §1a](../../spec/01-pipeline-flow.md#1a-engineering-constraints),
each row targets exactly one organelle via `assembly_target ∈
{animal_mt, plant_pt, plant_mt}`; a plant sample requiring both
plastid and mitogenome assemblies submits two rows sharing the same
reads:

- `TEST-PLANT-01-pt` — `assembly_target=plant_pt`, uses `plant.fastq.gz`
- `TEST-PLANT-01-mt` — `assembly_target=plant_mt`, uses `plant.fastq.gz`
- `TEST-ANIMAL-01` — `assembly_target=animal_mt`, uses `animal.fastq.gz`
- `TEST-ANIMAL-02` — `assembly_target=animal_mt`, uses
  `animal.fastq.gz|animal_b.fastq.gz` (exercises PARSE_SAMPLESHEET's
  multi-file concat path)

[`samples_bad.csv`](samples_bad.csv) — negative-path fixture with a
duplicate `sample_id` for `-profile test_bad_samplesheet`.

## Tier 2 (integration) fixtures

Not present yet. Planned: larger subsets (~5–50 MB) with enough coverage
to exercise real assembly, hosted on Azure blob and fetched by a
`fetch_integration_fixtures.sh` script. See the plan file for details.
