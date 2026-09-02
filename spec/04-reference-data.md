## 4. Reference data

Target-keyed reference bundle, versioned and configurable via
`params.organelle_refs` (bundle root; per-target `.mmi` indices live
inside):

| Assembly target | Organelle ref (recruit/bin) | Protein panel (miniprot) | Barcode selector |
|---|---|---|---|
| `plant_pt` | GetOrganelleDB `embplant_pt` (plastid) | ~79 plastid PCGs | rbcL, matK, ndhF, atpB, psbA, rpoB |
| `plant_mt` | RefSeq Viridiplantae mitochondrion, plastid-masked ([§4.1a](#41a-plant_mt-a-plastid-masked-refseq-panel)) | ~40 plant mt PCGs | cox1, cob, nad1, atp1, matR |
| `animal_mt` | GetOrganelleDB `animal_mt` (metazoan mitogenome) | 13 metazoan mt PCGs | COX1, CYTB, COX2, COX3, ND1, ATP6 |

The last two columns are **different things** and [§4.3](#43-protein-panel-and-barcode-selector)
explains why keeping them apart matters: the panel is what miniprot aligns,
the selector is which of those genes count as barcodes. Selector content is
**not authoritative here** — it is parsed at runtime from
[`assets/loci.json`](../assets/loci.json) ([brief.md §3.7](brief.md)); the
table shows representative loci only.

### 4.1 Setup task: source reference panel from GetOrganelle

Decision: bootstrap the per-target organelle references from
[GetOrganelleDB](https://github.com/Kinggerm/GetOrganelleDB) (the reference
library shipped with the GetOrganelle assembler). The libraries are already
target-partitioned and curated for organelle recruitment. We reuse them as
recruitment references without using GetOrganelle as the assembler.

**Mapping from GetOrganelle libraries to our assembly targets:**

| Assembly target | GetOrganelle library | Notes |
|---|---|---|
| `plant_pt` | `embplant_pt` | Plastid. Use `other_pt` in addition if non-embryophyte plant lineages need coverage. |
| `plant_mt` | ~~`embplant_mt`~~ | **Superseded — see [§4.1a](#41a-plant_mt-a-plastid-masked-refseq-panel).** GetOrganelleDB ships exactly one plant-mitochondrion seed, which is too thin for the three jobs we ask a panel to do. |
| `animal_mt` | `animal_mt` | Metazoan mitogenome. |

**Setup procedure (one-off, scripted, versioned):**

1. Fetch a pinned [GetOrganelleDB release](https://github.com/Kinggerm/GetOrganelleDB/releases) tarball; record the release tag in `params.organelle_refs.version`.
2. Extract per-library FASTA files.
3. Copy each library FASTA to `<target>.fa` — **no concatenation across targets** ([spec §1a](01-pipeline-flow.md#1a-engineering-constraints)); each `assembly_target` recruits against its own index.
4. Build minimap2 indices: `minimap2 -d <target>.mmi <target>.fa` for each of `plant_pt`, `plant_mt`, `animal_mt`.
5. Stage indices + source manifests under a versioned directory (e.g. `refs/v2026.06/recruit/`); point `params.organelle_refs` at that directory.
6. Record source URLs, release tags, and SHA256 digests in a `manifest.json` alongside the indices. Emit the manifest version into per-run metadata output ([brief.md §5](brief.md)).

**Re-evaluated after P1 — resolved for `plant_mt`, held for the other two.**
The clause above asked whether GetOrganelle library coverage would prove
insufficient. For `plant_mt` it did, decisively, and `plant_mt` now uses the
RefSeq-derived fallback this clause anticipated (§4.1a). `plant_pt` (101
genomes) and `animal_mt` (1853 genomes) remain GetOrganelle seed copies and
have shown no comparable symptom.

### 4.1a `plant_mt`: a plastid-masked RefSeq panel

The `plant_mt` panel is not a GetOrganelle seed copy. GetOrganelleDB ships a
**single** plant-mitochondrion seed (one *Vigna radiata* mitogenome, 398 kb).
GetOrganelle can afford that because its assembler grows outward from the seed
iteratively; we use the panel as a one-shot filter and — since the coverage
split at [§2.1.5](02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)
and the sibling discrimination at
[§3.7.1](03-organelles.md#371-homology-is-measured-against-sibling-panels-too)
— as the reference two decisions are made against. For those jobs one distant
genome is far too thin: genuine *Datura* mitochondrial contigs matched the
*Vigna* seed over only 20–34 % of their length.

**Why not the raw RefSeq set.** Plant mitogenomes contain NUPTs — plastid
sequence copied into the mitochondrial genome over evolutionary time. Across
RefSeq Viridiplantae these are present in most genomes, so a broad,
*unmodified* plant-mitochondrion panel is in part a chloroplast panel: plastid
contigs match the mitochondrial reference as well as they match the chloroplast
one, and the signal that separates the two organelles collapses. Measured on
`INT-PLANT-01-mt`, raw RefSeq drove the plastid sibling margin to ±0.000 and
re-opened a false `ok` at the coverage gate.

**The derivation.** Take the RefSeq Viridiplantae mitochondrion set, align it
against the `plant_pt` panel, merge the aligned query intervals per genome
([§3.7.2](03-organelles.md#372-aligned-fraction-is-merged-across-all-alignment-blocks)),
and replace them with `N`. Masking substitutes, never deletes, so coordinates
are preserved. Genomes measuring majority-plastid by this rule are dropped
rather than emitted as mostly-`N` references.

Masking operates on the **reference**, not on reads, so it is not a
[principle 5](../CONSTITUTION.md) concern: no read is excluded for resembling
anything. We are removing sequence from the bait that was never
mitochondrial-specific, so that "aligns to the plant mitochondrion panel"
means what it claims.

Implemented by [`scripts/refdata/mask_panel.py`](../scripts/refdata/mask_panel.py),
driven from `build_recruit.sh`. The masking parameters (alignment preset,
minimum alignment length, maximum masked fraction per genome) are script
arguments, not constants, and the values actually applied are recorded in
`manifest.json` along with the masked fraction and every skipped genome — see
task 28_plastid_masked_mt_panel.md for the measured evidence and the parameter
rationale.

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

### 4.3 Protein panel and barcode selector

**Two distinct things live here, and conflating them is the mistake this
section exists to prevent:**

1. **The protein panel** (`proteins/<origin>/<gene>.faa` in the reference
   bundle) — the sequences miniprot aligns. Since [§8 item 3](07-open-questions.md#8-remaining-open-questions)
   it is the organelle's **full protein-coding complement**: 13 genes for
   `animal_mt`, ~79 for `plant_pt`, ~40 for `plant_mt`. `MINIPROT_CDS`
   (stage 12) concatenates the origin's files into one query and produces
   `cds.gff`.
2. **The barcode selector** (`assets/loci.json`, `params.locus_panel`) —
   which of those genes are *barcodes*. `EXTRACT_BARCODES` (stage 13)
   picks the `cds.gff` features whose gene symbol appears here.

The panel is reference data; the selector is versioned config coupled to
Taxodactyl's accepted-loci set ([principle 9](../CONSTITUTION.md)). They
change for different reasons and on different cadences. The selector's
content and schema are **unchanged** by the move to a comprehensive panel
— it is now purely a statement about which loci are barcodes, rather than
doubling as a statement about which proteins get aligned.

**Consequence for the panel build:** a gene may be in the panel with no
usable representative sequence, and the completeness yardstick in
`assets/organelle_gene_sets.json` should still count it missing. The
yardstick is a third, separate list — expected gene *names* — and is not
the panel.

We maintain our own wf5-specific locus selector rather than inheriting
Taxodactyl's
config: Taxodactyl's `loci.json` is designed for **matching gene names
against annotation output** (via ambiguous/non-ambiguous synonyms) and mixes
bacterial, fungal, and nuclear loci that are out of scope for organelle
recovery. Our panel groups loci by **organelle origin** and stores only the
canonical NCBI gene symbol — the shape the barcode selector and the RefSeq
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

**Fetch behaviour (task 29):**

- **Source of protein sequences:** CDS translations parsed directly out of the RefSeq organelle GenBank flat files (`plastid.*.genomic.gbff.gz`, `mitochondrion.1.genomic.gbff.gz`) — the same RefSeq release already downloaded for [§4.2](#42-refseq-organelle-dbs-for-blast_validate)'s genomic FASTA, not fresh per-gene NCBI queries. Fewer moving parts (one bulk download + parse instead of ~130 rate-limited esearch/efetch round trips) and the provenance is already recorded. Kingdom split for the mitochondrion GBFF (not pre-split, unlike plastid) uses the GenBank record's own taxonomy lineage (`Metazoa` / `Viridiplantae`), not a fresh accession→taxid lookup. Parsing runs inside the `neoformit/daff-wf5-scripts:test` image (biopython pinned) via [`scripts/refdata/parse_gbff_cds.py`](../scripts/refdata/parse_gbff_cds.py); building the per-gene panel from the parsed CDS records is host-side, stdlib-only, in [`scripts/refdata/build_proteins.py`](../scripts/refdata/build_proteins.py) (`full` mode).
- **Gene list:** the organelle's full protein-coding complement, not the barcode subset (see the two-things note above) — the canonical list per target lives in [`assets/organelle_gene_sets.json`](../assets/organelle_gene_sets.json) (13 / ~79 / ~40 genes for `animal_mt` / `plant_pt` / `plant_mt`), with a `protein_coding_aliases` table mapping the canonical symbol to the `/gene` qualifier spellings actually seen in RefSeq (case-insensitive; an `MT-` locus prefix, as used by vertebrate mitogenome records, is stripped before matching). `plant_mt` is built from **mitochondrial CDS records only** — plastid-derived proteins leaking in would let a NUPT insertion be annotated as a functional mitochondrial gene, the same class of error [§4.1a](#41a-plant_mt-a-plastid-masked-refseq-panel) masks out of the recruitment panel. Do not mask the protein *sequences* (a mitochondrial protein is a mitochondrial protein); the discipline is entirely in which source records are used.
- **Per-gene curation:** 5–10 representatives per gene, chosen to maximise **genus spread** among the candidate CDS translations (a phylogenetic-diversity proxy, not a divergence measure) rather than to maximise count — miniprot tolerates divergence, and the marginal sensitivity gain from a 50th representative is small next to the first five spanning the clade ([§9 item 11](07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)). A hard minimum (≥3 candidates) excludes a gene from the panel if RefSeq annotation coverage is too thin; the gene still appears in `organelle_gene_sets.json`'s completeness yardstick, just absent from the `.faa` panel.
- **Version pinning:** record per-gene candidate/representative counts, chosen accessions, RefSeq release, and the parse method in `refs/<version>/proteins/provenance.json`, folded into `manifest.json` by `build_manifest.py`.
- **Build:** one FASTA per gene per origin: `proteins/<origin>/<gene>.faa`. Per-gene files keep the barcode subset trivially selectable and the bundle diff-able; `MINIPROT_CDS` concatenates them at runtime.
- **Panel-competition check.** Broadening the query set changes what competes at each locus — with ~790 proteins in play instead of ~60, miniprot may select a different best alignment at a barcode locus. Since barcodes are the contractual output, any change to panel breadth must be measured before it ships: run miniprot over the same target contigs with the old and new panels and compare, per barcode locus, whether the locus is still called, whether coordinates moved, and whether identity and query coverage improved. **No barcode locus lost, none degraded** — a coordinate shift is acceptable only where the metrics improve. Run and recorded for `v2026.09` in task 29's Outcomes.
- **Refresh trigger:** re-fetch when the gene list or the barcode selector changes, or annually, whichever comes first.

### 4.4 Consolidated build script

All three reference-build tasks (§4.1–§4.3) live in a single versioned script
(`scripts/build_refs.sh` or equivalent) that emits `refs/v<YYYY.MM>/` with:

```
refs/v2026.08/
├── recruit/
│   ├── plant_pt.fa       # GetOrganelle embplant_pt (plastid)
│   ├── plant_pt.mmi
│   ├── plant_mt.fa       # RefSeq Viridiplantae mt, plastid-masked (§4.1a)
│   ├── plant_mt.mmi
│   ├── plant_mt.masking.json   # masking params + per-genome fractions
│   ├── animal_mt.fa      # GetOrganelle animal_mt
│   └── animal_mt.mmi
├── validate/
│   ├── refseq_pt.{nhr,nin,nsq,...}      # plant plastid BLAST DB
│   ├── refseq_mt_metazoa.{...}
│   ├── refseq_mt_viridiplantae.{...}
│   └── refseq_mt_viridiplantae.fa       # retained: input to §4.1a masking
├── proteins/             # full protein-coding complement per origin (§4.3)
│   ├── animal_mt/<GENE>.faa      # 13 genes
│   ├── plant_pt/<gene>.faa       # ~79 genes
│   └── plant_mt/<gene>.faa       # ~40 genes, mitochondrial CDS only
├── annotate/
│   └── mitos/refseq89m/  # MITOS2 reference data (metazoan, RefSeq 89)
└── manifest.json         # SHA256 + version + source URL for every input
```

`manifest.json` is the single source of truth for reference provenance; its
version string is emitted into per-run metadata ([brief.md §5](brief.md)) so
every result is traceable to the exact reference bundle used.

`params.recruit_panels` / `protein_panel` / `blast_db` / `annotate_refs`
each point *inside* `refs/<version>/` — at a subdirectory, not the bundle
root — so none of them can locate `manifest.json` itself.
`params.refs_manifest` is the bundle-root param that does: it points
directly at `refs/<version>/manifest.json` and is staged into `COLLATE`
(C6) as a Nextflow `path` input (pattern 1 of
[§1a](01-pipeline-flow.md#1a-staging-reference-files), required by
[rule 13](../CONSTITUTION.md) — a bare `${params.x}` string is not
stageable on a remote executor). `COLLATE` copies only `version` and
`generated_at` into each sample's `metadata.json`; the full manifest is
`RUN_REPORT`'s (C7) to inline into `run_manifest.json` (task 42 §5.2).

### 4.5 Fetching and verifying a bundle

`scripts/fetch_refs.sh <version>` downloads `refs-<version>.tar.gz` from Azure
blob, verifies the tarball's SHA256 against `scripts/refs-<version>.sha256`,
and unpacks it into `refs/<version>/`. The SHA256 check only proves the
*tarball* arrived intact — it says nothing about whether the tarball was
*built* with every component `manifest.json` lists inside it (task 41's
incident: a bundle built before `annotate/` existed passed its checksum
cleanly and failed 12 minutes later, inside a container, on a dangling
`annotate/` symlink). So after unpacking, `fetch_refs.sh` also reads the
unpacked `manifest.json`'s `artefacts` keys, derives the set of top-level
components they belong to (`recruit`, `validate`, `proteins`, `annotate`,
...), and fails with the missing names if any of them is not present on
disk. This is provenance-tracking applied one layer up from the manifest
itself ([principle 10](../CONSTITUTION.md)) — a missing input must never be
reported as a successful fetch.
