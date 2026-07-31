# Task 2 — Reference data bundle

**Phase:** P1 prerequisite (feeds real logic in [plan.md §6](../../plan.md)).
**Goal:** A consolidated, versioned reference bundle `refs/v<YYYY.MM>/`
containing every offline reference dataset consumed by the pipeline —
`recruit/` (§4.1), `validate/` (§4.2), `proteins/` (§4.3) — with a
`manifest.json` that pins provenance for each input.

**Exit criteria:**

- `scripts/build_refs.sh` (or `.py`) produces `refs/v2026.07/` matching the
  tree in [plan.md §4.4](../../spec/04-reference-data.md#44-consolidated-build-script).
- Every artefact recorded in `manifest.json` with source URL, upstream
  release/version tag, download date, and SHA256 digest.
- The bundle is a drop-in target for `params.kingdom_refs`,
  `params.blast_db`, `params.locus_panel` — the P0 stub workflow points at
  it and stages files correctly under `-profile test`.
- Build runs in a single container image (`wf5/ref-build:<tag>`) — no host
  tooling assumed.

**Not in scope:**

- Uploading the bundle to Azure blob storage (that belongs to
  [`deploy/azure/`](../../deploy/azure/), triggered by a bundle-version bump —
  see [docs/azure/README.md](../../docs/azure/README.md)).
- Any pipeline-side change; this task ends when the bundle exists on disk
  and validates against the manifest.

**Cross-cutting rules:**

- **Reproducible from source.** No manual edits to fetched files. Every
  transformation (concat, kingdom-split, index build) is scripted.
- **Pin everything.** Upstream release tag or dated snapshot for every
  input; SHA256 on the final artefact. A rebuild from the same script on
  the same day should produce byte-identical output modulo timestamps.
- **Single container.** The build image bakes `minimap2`, `blast+`,
  `seqkit`, `entrez-direct` / `datasets`, `pandas`, `curl`,
  `taxonkit`. Match the "deps in image, code at runtime" pattern
  ([plan.md §2.2](../../plan.md)) — the build script lives under `scripts/`
  and is mounted, not baked.

---

## 1. Target layout

```
refs/v2026.07/
├── recruit/
│   ├── plant_pt.fa                    # embplant_pt (plastid)
│   ├── plant_pt.mmi                   # minimap2 -d
│   ├── plant_mt.fa                    # embplant_mt (mitogenome)
│   ├── plant_mt.mmi
│   ├── animal_mt.fa                   # animal_mt
│   └── animal_mt.mmi
├── validate/
│   ├── refseq_pt.{nhr,nin,nsq,...}    # plastid BLAST DB (plants)
│   ├── refseq_mt_metazoa.{...}        # split from refseq mitochondrion
│   └── refseq_mt_viridiplantae.{...}
├── proteins/
│   ├── animal_mt/<GENE>.faa           # canonical NCBI gene symbols (see assets/loci.json)
│   ├── plant_pt/<gene>.faa
│   └── plant_mt/<gene>.faa
└── manifest.json
```

`scripts/build_refs.sh` writes to a staging dir, verifies checksums, then
atomically renames into place. The bundle root is immutable once written —
a new version bump means a new `v<YYYY.MM>/` directory.

### 1.1 Disk footprint (estimate)

Final bundle — what gets published to blob storage and staged onto compute
nodes:

| Path | Estimated size |
|---|---|
| `recruit/` (FASTAs + .mmi indices) | ~100 MB |
| `validate/refseq_pt.*` (BLAST DB) | ~10 GB |
| `validate/refseq_mt_metazoa.*` | ~4 GB |
| `validate/refseq_mt_viridiplantae.*` | ~1 GB |
| `proteins/` (per-locus FASTAs) | <100 MB |
| **Bundle total** | **~15 GB** |

Working set during a build (raw downloads + intermediates, cleaned up
before the atomic rename) roughly doubles that:

| Item | Estimated size |
|---|---|
| RefSeq plastid gz + unpacked FASTAs | ~8 GB |
| RefSeq mitochondrion gz + unpacked | ~4 GB |
| `nucl_gb.accession2taxid.gz` (extracted) | ~10 GB |
| taxdump (already staged at `/home/cameron/.taxonkit/`) | ~500 MB |
| Kingdom-split intermediate FASTAs before `makeblastdb` | ~5 GB |
| **Peak scratch** | **~30 GB** |

Provision **~50 GB free** on the build host to be safe. Numbers scale
roughly linearly with RefSeq release growth (~5% per quarter); re-check
before a bundle bump.

## 2. Inputs & sources

| Component | Source | Version pin |
|---|---|---|
| `recruit/plant_pt.fa` | [GetOrganelleDB](https://github.com/Kinggerm/GetOrganelleDB) release `0.0.1`, `SeedDatabase/embplant_pt.fasta` | GitHub release tag |
| `recruit/plant_mt.fa` | GetOrganelleDB release `0.0.1`, `SeedDatabase/embplant_mt.fasta` | GitHub release tag |
| `recruit/animal_mt.fa` | GetOrganelleDB release `0.0.1`, `SeedDatabase/animal_mt.fasta` | GitHub release tag |
| `validate/refseq_pt.*` | [`ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*.genomic.fna.gz`](https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/) | RefSeq release number (`RELEASE_NUMBER` file) + download date |
| `validate/refseq_mt_*` | [`ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*.genomic.fna.gz`](https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/), kingdom-split via NCBI taxonomy | RefSeq release + taxdump date |
| `proteins/<origin>/<gene>.faa` | [`assets/loci.json`](../../assets/loci.json) (origin-keyed panel of canonical NCBI gene symbols) → per-locus RefSeq protein queries | `loci.json` git SHA + query date |

### 2.1 GetOrganelleDB (recruit)

Already downloaded locally under
`reference-material/getorganelle-db/0.0.1/`. The build script accepts a
local path (`--getorganelle-db <path>`) to avoid re-downloading, but must
also support fetching if absent. Use `SeedDatabase/` — that is what
GetOrganelle uses to seed recruitment; `LabelDatabase/` is downstream
labelling and not needed here.

Steps — one FASTA + index per `<kingdom, organelle>` pair; **no
concatenation** ([plan.md §4.1](../../spec/04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle)):
plant cp and plant mt stay as distinct panels so downstream stages can
distinguish per-contig which organelle a hit came from.

1. Copy `embplant_pt.fasta` → `recruit/plant_pt.fa`.
2. Copy `embplant_mt.fasta` → `recruit/plant_mt.fa`.
3. Copy `animal_mt.fasta` → `recruit/animal_mt.fa`.
4. Build indices: `minimap2 -d recruit/<name>.mmi recruit/<name>.fa` for
   each of `plant_pt`, `plant_mt`, `animal_mt`.

### 2.2 RefSeq organelle (validate)

Two RefSeq collections; mitochondrion needs a kingdom split
([plan.md §4.2](../../plan.md)):

1. Mirror `refseq/release/plastid/*.genomic.fna.gz` and
   `refseq/release/mitochondrion/*.genomic.fna.gz`; note release number
   from `RELEASE_NUMBER`.
2. Fetch `taxdump.tar.gz` (already downloaded to
   /home/cameron/.taxonkit/new_v0.20/) + `nucl_gb.accession2taxid.gz`; note
   dump date.
3. For mitochondrion: parse each FASTA header's accession, look up taxid,
   walk taxonomy up to superkingdom (you can use
   /home/cameron/.local/bin/taxonkit for ths) → route to
   `refseq_mt_metazoa.fa` or `refseq_mt_viridiplantae.fa`. Records
   outside those two kingdoms are logged and discarded (fungi, protists,
   etc — not in scope for this pipeline).
4. `makeblastdb -in refseq_<...>.fa -dbtype nucl -parse_seqids -out
   validate/refseq_<...>` for all three DBs.

### 2.3 Protein panel (proteins)

Locus list is authoritative in [`assets/loci.json`](../../assets/loci.json)
(schema `wf5/loci-panel/v1`) — a wf5-specific panel grouped by organelle
origin, holding only canonical NCBI gene symbols. We do **not** inherit
Taxodactyl's `loci.json` here: Taxodactyl's schema is for matching gene
names against annotation output (via ambiguous/non-ambiguous synonyms) and
mixes bacterial/fungal/nuclear loci that are out of scope
([plan.md §4.3](../../spec/04-reference-data.md#43-protein-panel-for-miniprot_extract)).

Panel shape:

```json
{
  "animal_mt": ["COX1", "COX2", "COX3", "CYTB", "ND1", "ATP6"],
  "plant_pt":  ["rbcL", "matK", "atpB", "ndhF", "psbA", "rpoB"],
  "plant_mt":  ["cox1", "cob", "nad1", "atp1", "matR"]
}
```

The build script:

1. Reads the panel from `--loci-config <path>` (same file
   `params.locus_panel` points at — single source of truth).
2. For each `<origin, gene>` entry, queries NCBI protein via
   `esearch → efetch` (POST, to avoid HTTP 414 on long ID batches).
   Origin → taxid restriction:
   - `animal_mt` → `txid33208[Organism]` (Metazoa)
   - `plant_pt` / `plant_mt` → `txid33090[Organism]` (Viridiplantae)
3. Writes `proteins/<origin>/<gene>.faa`, capped at
   `--max-per-locus` (default 50).
4. Enforces a per-locus minimum (default 3) — fails the build if any
   locus returns too few hits, which would silently break
   `MINIPROT_EXTRACT`.

Query strategy per locus: `<GENE>[Gene Name] AND <origin-taxid>[Organism]
AND refseq[Filter]`. NCBI's `[Gene Name]` field is case-insensitive but
expects **gene symbols** (e.g. `COX1`), not descriptions
(`"cytochrome oxidase subunit 1"` returns zero hits). Record the exact
query string in the manifest entry so re-running is deterministic.

## 3. `manifest.json` schema

```json
{
  "bundle_version": "v2026.07",
  "built_at": "2026-07-26T14:30:00Z",
  "built_by": "cameron@neoformit.com",
  "builder_image": "wf5/ref-build:0.1.0",
  "artefacts": {
    "recruit/plant_pt.fa": {
      "sha256": "…",
      "sources": [
        {
          "url": "https://github.com/Kinggerm/GetOrganelleDB/releases/tag/0.0.1",
          "path_in_source": "SeedDatabase/embplant_pt.fasta",
          "release_tag": "0.0.1"
        }
      ]
    },
    "validate/refseq_mt_metazoa.nsq": {
      "sha256": "…",
      "sources": [
        {
          "url": "https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/",
          "refseq_release": "227",
          "downloaded_at": "2026-07-26",
          "taxdump_date": "2026-07-20"
        }
      ]
    },
    "proteins/animal_mt/COX1.faa": {
      "sha256": "…",
      "sources": [
        {
          "url": "https://www.ncbi.nlm.nih.gov/protein/",
          "query": "\"COX1\"[Gene Name] AND txid33208[Organism] AND refseq[Filter]",
          "downloaded_at": "2026-07-26",
          "record_count": 42
        }
      ]
    }
  }
}
```

The `bundle_version` field is what `RUN_REPORT` (C7) emits into
`run_manifest.json` for provenance ([plan.md §2 stage 16](../../plan.md)).

## 4. Script structure

`scripts/build_refs.sh` — orchestrator; delegates each sub-build to a
helper under `scripts/refdata/`:

```
scripts/
├── build_refs.sh                # orchestrator; parses args, calls helpers
└── refdata/
    ├── build_recruit.sh         # §2.1 GetOrganelleDB → recruit/
    ├── build_validate.sh        # §2.2 RefSeq → validate/
    ├── build_proteins.py        # §2.3 Taxodactyl config → proteins/
    ├── split_refseq_mt.py       # taxonomy-based kingdom split
    └── write_manifest.py        # walks bundle dir, computes SHA256s
```

Orchestrator flags:

```
--version v2026.07            # required; bundle version string
--out refs/                   # bundle root (default: ./refs)
--getorganelle-db PATH        # optional; skip download if present
--taxodactyl-config PATH      # required for proteins/
--only recruit|validate|proteins  # optional; skip sub-builds for iteration
--dry-run                     # print planned actions, fetch nothing
```

Fail fast on any sub-build error; a partial bundle is never renamed into
`refs/v<...>/`.

## 5. Builder container

Image: `wf5/ref-build:<tag>` under `containers/ref-build/Dockerfile`.
Base: `python:3.12-slim`. Installs:

- `minimap2` (apt)
- `ncbi-blast+` (apt)
- `seqkit` (github release)
- `entrez-direct` **or** `ncbi-datasets-cli` — pick one in §7
- `taxonkit` (github release, needs taxdump path)
- Python: `pandas`, `biopython`, `requests`

Image tag bumps on Dockerfile / requirements change only, per
[plan.md §5b](../../plan.md).

## 6. Deliverables checklist

- [ ] `scripts/build_refs.sh` + helpers under `scripts/refdata/`
- [ ] `containers/ref-build/Dockerfile` + built image published to
      `neoformit/daff-wf5-ref-build:<tag>`
- [ ] `refs/v2026.07/` on the build machine (**not** committed to git —
      too big; add `refs/` to `.gitignore`)
- [ ] `manifest.json` validates against §3 schema
- [ ] P0 stub workflow runs successfully against
      `--kingdom_refs refs/v2026.07/recruit`,
      `--blast_db refs/v2026.07/validate`,
      `--locus_panel refs/v2026.07/proteins` (i.e. the split-directory
      layout that `run-wf5.sh` already expects)
- [ ] `docs/reference-data.md` updated with build invocation + refresh
      cadence

## 7. Open questions

1. **`entrez-direct` vs `datasets`** for protein fetching. Both work;
   `datasets` has cleaner JSON output and is Anthropic-friendly for
   scripting, `entrez-direct` is battle-tested and closer to the query
   strings we already write. Recommend **`datasets`** unless a locus
   query is easier to express as an Entrez `esearch | efetch` pipe.
2. **Where to run the build.** Local workstation vs an Azure Batch job
   pointed at the same pool. Recommend local for the first bundle
   (fast iteration, no batch queue latency) and revisit if the RefSeq
   downloads become painful over the office link.
3. **Refresh cadence and version scheme.** Plan §4.2 says quarterly for
   RefSeq. Suggest bundle version = year-month of the RefSeq release
   used (`v2026.07` for RefSeq 227), with a bump whenever any input
   changes — a fresh Taxodactyl locus config mid-quarter still triggers
   a new bundle even if RefSeq hasn't rolled.
4. **Committing the manifest.** `refs/` is git-ignored, but the
   `manifest.json` for each released bundle is small and worth keeping
   for audit — propose committing `refs/manifests/v<YYYY.MM>.json` even
   though the artefacts themselves live in blob storage.
