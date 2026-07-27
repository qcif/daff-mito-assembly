## 4. Reference data

Kingdom-keyed reference bundle, versioned and configurable via `params.kingdom_refs`:

| Kingdom | Organelle refs | Protein panel (for miniprot) |
|---------|---------------|------------------------------|
| plant | RefSeq plastid + plant mitogenome | rbcL, matK, ndhF, atpB, COX1, etc. |
| animal | RefSeq metazoan mitogenome | COX1, CYTB, 12S, 16S, etc. |

Locus panel content is **not authoritative here** — it is parsed at runtime
from Taxodactyl's accepted-loci config ([brief.md §3.7](brief.md)). The table
above shows representative loci only.

### 4.1 Setup task: source reference panel from GetOrganelle

Decision: bootstrap the kingdom organelle reference panel from
[GetOrganelleDB](https://github.com/Kinggerm/GetOrganelleDB) (the reference
library shipped with the GetOrganelle assembler). The libraries are already
kingdom-partitioned and curated for organelle recruitment. We reuse them as
recruitment references without using GetOrganelle as the assembler.

**Mapping from GetOrganelle libraries to our kingdom panels:**

| Our kingdom | GetOrganelle libraries | Notes |
|---|---|---|
| plant | `embplant_pt` + `embplant_mt` | Plastid + mitogenome. Use `other_pt` if non-embryophyte plant lineages need coverage. |
| animal | `animal_mt` | Mitogenome only. |

**Setup procedure (one-off, scripted, versioned):**

1. Fetch a pinned [GetOrganelleDB release](https://github.com/Kinggerm/GetOrganelleDB/releases) tarball; record the release tag in `params.kingdom_refs.version`.
2. Extract per-library FASTA files.
3. Concatenate the libraries listed above into one FASTA per kingdom.
4. Build minimap2 indices: `minimap2 -d kingdom_<plant|animal>.mmi kingdom_<...>.fa`.
5. Stage indices + source manifests under a versioned directory (e.g. `refs/v2026.06/`); point `params.kingdom_refs` at it.
6. Record source URLs, release tags, and SHA256 digests in a `manifest.json` alongside the indices. Emit the manifest version into per-run metadata output ([brief.md §5](brief.md)).

**Re-evaluate after P1.** If GetOrganelle library coverage proves insufficient on the P1 plant test data (older or non-model lineages underrepresented), fall back to a self-built RefSeq-derived panel using the canonical sources listed in §4.

### 4.2 RefSeq organelle DBs (for `BLAST_VALIDATE`)

Stage 11 needs kingdom-appropriate BLAST nucleotide DBs, separate from the
recruitment panel — recruitment is fuzzy and permissive; validation must be
against a clean, curated reference set to catch mis-binning.

- **Source:** NCBI RefSeq organelle FTP:
  - plastid: `https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*.genomic.fna.gz`
  - mitochondrion: `https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*.genomic.fna.gz`
- **Kingdom split:** RefSeq mitochondrion is not pre-split by kingdom. Split by NCBI taxonomy at build time using accession → taxid mapping (`efetch` or the `nucl_gb.accession2taxid` dump) and the taxonomy dump (`taxdump.tar.gz`). Emit two files: `refseq_mt_metazoa.fa`, `refseq_mt_viridiplantae.fa`. Plastid is plants-only by definition.
- **Version pinning:** record RefSeq release number (e.g. `RefSeq 227`) + download date + taxdump date in `manifest.json`.
- **Build:** `makeblastdb -in refseq_<kingdom>_<organelle>.fa -dbtype nucl -parse_seqids -out blastdb/refseq_<...>`.
- **Refresh cadence:** re-fetch quarterly; RefSeq updates every ~2 months.

### 4.3 Protein panel for `MINIPROT_EXTRACT`

miniprot needs a protein FASTA per locus (protein-to-genome alignment). We
maintain our own wf5-specific locus panel rather than inheriting Taxodactyl's
config: Taxodactyl's `loci.json` is designed for **matching gene names
against annotation output** (via ambiguous/non-ambiguous synonyms) and mixes
bacterial, fungal, and nuclear loci that are out of scope for organelle
recovery. Our panel groups loci by **organelle origin** and stores only the
canonical NCBI gene symbol — the shape MINIPROT_EXTRACT and the RefSeq
protein query need.

**Schema** (`assets/loci.json`, referenced by `params.locus_panel`):

```json
{
  "$schema": "wf5/loci-panel/v1",
  "animal_mt": ["COX1", "COX2", "COX3", "CYTB", "ND1", "ATP6"],
  "plant_pt":  ["rbcL", "matK", "atpB", "ndhF", "psbA", "rpoB"],
  "plant_mt":  ["cox1", "cob", "nad1", "atp1", "matR"]
}
```

- **Origin keys** match the [§4.1](#41-setup-task-source-reference-panel-from-getorganelle) recruit-panel naming (`plant_pt`, `plant_mt`, `animal_mt`) so downstream stages route consistently.
- **Values** are canonical NCBI `[Gene Name]` symbols. NCBI is case-insensitive on this field.
- **Protein-coding only.** rRNA barcodes (12S, 16S) and tRNA barcodes (trnL) are excluded because miniprot needs protein; those would need a separate nucleotide-based extraction path (not scoped in P1).
- The file is staged onto every consuming process via a value channel and parsed inside the container ([§1a](01-pipeline-flow.md#1a-engineering-constraints) channel rule); it is never read from `params.*` inside a `script:` block or preloaded into a Groovy map at workflow entry.

**Fetch behaviour:**

- **Source of protein sequences:** [NCBI RefSeq protein](https://ftp.ncbi.nlm.nih.gov/refseq/release/) — query by gene symbol + kingdom taxid restriction. Kingdom mapping is derived from the origin key: `animal_mt → txid33208[Organism]` (Metazoa), `plant_pt` / `plant_mt → txid33090[Organism]` (Viridiplantae). Alternative: [UniProt](https://www.uniprot.org/) reviewed entries.
- **Per-locus curation:** pull up to ~50 representative proteins per locus; miniprot tolerates divergence — 5–10 representatives per locus is enough; more slows alignment without improving sensitivity. A hard minimum (≥3) fails the build if a locus returns too few hits.
- **Version pinning:** record per-locus accessions + fetch date + exact query string in `manifest.json` under `locus_panel_proteins`.
- **Build:** one FASTA per locus per origin: `proteins/<origin>/<locus>.faa`. Concat per origin for a single miniprot invocation, or pass individually depending on runtime characteristics measured in P3.
- **Refresh trigger:** re-fetch when the locus panel changes, or annually, whichever comes first.

### 4.4 Consolidated build script

All three reference-build tasks (§4.1–§4.3) live in a single versioned script
(`scripts/build_refs.sh` or equivalent) that emits `refs/v<YYYY.MM>/` with:

```
refs/v2026.06/
├── recruit/
│   ├── plant_pt.fa       # GetOrganelle embplant_pt (plastid)
│   ├── plant_pt.mmi
│   ├── plant_mt.fa       # GetOrganelle embplant_mt (mitogenome)
│   ├── plant_mt.mmi
│   ├── animal_mt.fa      # GetOrganelle animal_mt
│   └── animal_mt.mmi
├── validate/
│   ├── refseq_pt.{nhr,nin,nsq,...}      # plant plastid BLAST DB
│   ├── refseq_mt_metazoa.{...}
│   └── refseq_mt_viridiplantae.{...}
├── proteins/
│   ├── animal_mt/<GENE>.faa
│   ├── plant_pt/<gene>.faa
│   └── plant_mt/<gene>.faa
└── manifest.json         # SHA256 + version + source URL for every input
```

`manifest.json` is the single source of truth for reference provenance; its
version string is emitted into per-run metadata ([brief.md §5](brief.md)) so
every result is traceable to the exact reference bundle used.
