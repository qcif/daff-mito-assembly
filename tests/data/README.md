# Test fixtures

Real ONT reads for CI wiring tests. Small enough to commit; not large
enough to complete an actual organelle assembly — that's a Tier 2
integration concern, tracked separately.

## Files

| File | SHA256 | Reads | Size | Purpose |
|---|---|---|---|---|
| `animal.fastq.gz` | `fdc34cf1ae5bf8bdf4aa5cab08b4b5d915988d0a1f92df5a313954e5746ba000` | 70 | 502 K | `TEST-ANIMAL-01` primary; `TEST-ANIMAL-02` first file |
| `animal_b.fastq.gz` | `900c14fddd1b63efad4e6e487ea991e5aeed73de08c06daf358049047f476778` | 35 | 221 K | `TEST-ANIMAL-02` second file — different seed so read names don't collide with `animal.fastq.gz` (filtlong rejects duplicate read names) |
| `plant.fastq.gz` | `1e8d706a77e3d1cb66bd656ab8c26aa4115f211d258aace8a1b7076eba1e51c2` | 30 | 589 K | `TEST-PLANT-01-pt` and `TEST-PLANT-01-mt` (same reads, two `assembly_target` rows). Fewer reads than the animal fixture because Datura ONT reads are ~2× longer (~15 kb mean vs ~8 kb). |

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

```bash
ANIMAL_SRC="$PWD/tests/data/staging/SRR8306868_1.fastq.gz"
PLANT_SRC="$PWD/tests/data/staging/SRR11315861_1.fastq.gz"

# animal.fastq.gz — 70 reads, seed 42, ~500 KB target
docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v "$ANIMAL_SRC":/in.fastq.gz:ro \
    -v "$PWD/tests/data":/out \
    quay.io/biocontainers/seqtk:1.4--he4a0461_2 \
    bash -c 'seqtk sample -s42 /in.fastq.gz 70 | gzip > /out/animal.fastq.gz'

# animal_b.fastq.gz — 35 reads, seed 137 (non-overlapping with animal.fastq.gz)
docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v "$ANIMAL_SRC":/in.fastq.gz:ro \
    -v "$PWD/tests/data":/out \
    quay.io/biocontainers/seqtk:1.4--he4a0461_2 \
    bash -c 'seqtk sample -s137 /in.fastq.gz 35 | gzip > /out/animal_b.fastq.gz'

# plant.fastq.gz — 30 reads, seed 42, ~500 KB target
# (Datura reads are ~2× longer than aphid, so fewer reads to hit size.)
docker run --rm -u $(id -u):$(id -g) --platform linux/amd64 \
    -v "$PLANT_SRC":/in.fastq.gz:ro \
    -v "$PWD/tests/data":/out \
    quay.io/biocontainers/seqtk:1.4--he4a0461_2 \
    bash -c 'seqtk sample -s42 /in.fastq.gz 30 | gzip > /out/plant.fastq.gz'
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
