# miniprot — protein-to-genome alignment for `MINIPROT_EXTRACT`

**Reference:** Li H. (2023) *Protein-to-genome alignment with miniprot*.
Bioinformatics 39(1):btad014.
[doi:10.1093/bioinformatics/btad014](https://doi.org/10.1093/bioinformatics/btad014).
Source: [github.com/lh3/miniprot](https://github.com/lh3/miniprot).

## 1. What miniprot does

miniprot aligns **protein sequences to a complete genome**. It's the
protein-to-genome analogue of minimap2 (same author, same seed-chain-extend
strategy adapted for the protein/DNA case): index the genome, hash-map
k-mers over both strands in all six ORFs, chain matching anchors, then
close gaps with a dynamic-programming (DP) extension that models introns,
frameshifts, splice signals, and stop codons.

The novel bit vs. existing tools (GeneWise, Exonerate, Spaln2, GeMoMa):

- **Reduced amino acid alphabet** (SE-B(14), merging N/D) → compact
  4-bit encoding, dense hash index (~256/2 positions/entry, ~2^39 total
  for a human genome — fits in 32-bit ints).
- **Vectorised DP** using SIMD (SSE2 / ARM NEON), 8-way parallel with
  16-bit scores → **30–50× faster** than Spaln2/GeMoMa at comparable
  accuracy on real data (Table 1 of the paper).
- **Always aligns through introns** rather than skipping them — the DP
  is fast enough not to need the "localise-then-refine" trick older tools
  use. This is what gives miniprot its stability on distant homologs.

In numbers from the paper: mapping ~25k human proteins to the zebrafish
genome takes ~4 minutes (vs. ~2.5 minutes for Spaln2, ~10x faster than
GeMoMa). Base sensitivity 63.1% (Spaln2: 55.7%), specificity 94.9%.

## 2. Why it fits `MINIPROT_EXTRACT`

Stage 13 of the pipeline ([plan spec §2](04-reference-data.md#43-protein-panel-for-miniprot_extract),
[brief.md §7](../brief.md)) needs to locate a small panel of **canonical
barcode loci** (COX1, CYTB, rbcL, matK, …) on a *just-assembled* organelle
contig from a possibly-distant taxon.

Why protein-to-genome (not nucleotide-to-nucleotide) is the right shape:

| Requirement | Fit |
|---|---|
| Works across large evolutionary distance (no same-species reference) | ✓ — proteins are conserved where CDS nucleotides diverge. This is miniprot's headline use case (§4 Discussions: "Now with miniprot, we can perform approximate mapping and exact splice alignment in one go"). |
| Handles frameshifts in draft assemblies | ✓ — DP explicitly models frameshifts with a penalty (§2.6, Eq. 4–5). We can *use* frameshift calls as a QC signal. |
| Kingdom-appropriate genetic code | ✓ — `--trans` selects an NCBI translation table (2 = vertebrate mt, 5 = invertebrate mt, 4 = mould mt, 11 = plant plastid). Load-bearing for animal mt where table 2 vs. 5 changes stop/start codons. |
| Small reference input (a handful of protein sequences) | ✓ — the paper's tests are 20k+ proteins; we hand it ~50 proteins per locus. Runtime is trivial (seconds). |
| Precise CDS coordinates + protein alignment | ✓ — `--gff` emits GFF3 with CDS features; `--aln` includes the protein-to-CDS alignment. Both are needed for the barcode-extraction step. |

## 3. How miniprot works (brief)

1. **Index the genome** with `miniprot -d <ref.mpi> <ref.fa>` (or pass
   the FASTA directly to alignment; indexing is on-the-fly if omitted).
   Genome is scanned in all six ORFs; k-mers (default k = 6 aa) from
   translated ORFs ≥ 30 aa are hashed under the reduced alphabet.
2. **Seed** by scanning the query protein's k-mers against the index →
   list of (genome-pos, protein-pos) anchors.
3. **Chain** anchors by co-linear score, penalising gaps by an
   affine-with-log function that tolerates real introns
   (up to `--max-intron`, default 200 kb).
4. **Second chaining pass** over 5-mers around the top chains, at base
   resolution (no binning) — recovers small exons.
5. **DP extension** closes gaps between chained anchors and extends
   from terminal anchors. The DP has states for match/insertion/deletion
   *plus* three intron phases (A/B/C for phase-0/1/2 introns) and
   frameshift transitions (Eq. 7). Splice signals scored via a small
   position-weight model (default: GT-AG / GC-AG / AT-AC donor/acceptor
   with tuned penalties, §2.9).
6. **Pseudogene filter** re-scores single-exon alignments with a virtual
   intron cost, preferring the spliced hit when it's competitive.

## 4. Command shape for our use case

**Reference (genome):** the target contig from `METAFLYE` (Flye's
built-in polish already applied) — a single-record FASTA per sample.

**Query (protein panel):** the per-locus FASTA panel we build under
`refs/v<...>/proteins/<origin>/<gene>.faa`. Because the panel is small,
we can concat it per-organelle for one invocation per sample:

```bash
cat refs/v.../proteins/animal_mt/*.faa > proteins.animal_mt.faa

miniprot \
    -t 4 \
    --trans 2 \                # vertebrate mt genetic code (see §5.1)
    --gff \                    # GFF3 output
    --aln \                    # include protein-to-CDS alignment lines
    -I \                       # skip index build if genome tiny
    ${meta.sample_id}.assembly.fasta \
    proteins.animal_mt.faa \
    > ${meta.sample_id}.miniprot.gff
```

The GFF3 has one `mRNA` per hit with child `CDS` features (multi-exon if
plant mt has an intron; single-exon for most animal mt). `Target=<protein
id> <start> <end>` on the mRNA line carries the source protein accession —
we use that to route each hit back to its `<origin>/<gene>` and extract
the corresponding contig substring into `barcodes.fasta`.

## 5. Options that matter for wf5

### 5.1 `--trans <int>` — NCBI genetic code table

**Load-bearing** for animal mt. Options relevant to us:

| Table | Applies to | Notes |
|---|---|---|
| 1 | Standard (default) | plant nucleus — not applicable here |
| **2** | Vertebrate mt | UGA = Trp, AGR = stop |
| **5** | Invertebrate mt | UGA = Trp, AGR = Ser |
| 4 | Mould / protozoan mt | not scoped (fungal samples out of scope) |
| **11** | Plant plastid + bacterial | plant `plant_pt` panel |
| 1 | Plant mitochondrion | plant mt uses the *standard* code — no `--trans` override needed |

Handling in C5 (`validate_barcodes.py`, [plan §2.2](02-stages.md#22-custom-logic-components)):
kingdom-dispatch the table, and for animal mt run the **clade trial**
already planned in [§3.3](03-organelles.md#33-specific-issues-and-decisions):
try 2 then 5, pick the run whose ORF is intact, record the winning table
in `metadata.json`.

### 5.2 `--gff` vs `--gff-only` vs `--aln`

- `--gff` — emit GFF3 alongside the PAF-like default. **Use this.**
- `--gff-only` — GFF3 only, no PAF header.
- `--aln` — append alignment lines (`Prot:` / `Tran:` / `Geno:`) to the
  GFF. Needed for the ORF-integrity + internal-stop checks in C5.

### 5.3 `--outn <int>` — secondary alignments per query

Default is 100. For our short panel (~5–10 loci per organelle) we want
**one primary hit per query**, so `--outn 1` is safer — prevents multiple
overlapping hits per gene cluttering the extraction step. Secondaries
that miniprot considers viable would still surface as a QC warning in C5.

### 5.4 `--max-intron <int>`

Default 200 kb — irrelevant for animal mt (no introns) and mostly fine
for plant plastid. Plant mt can have long introns; leaving at default is
safe. Not a param we expect to tune.

### 5.5 `-I` / `-d`

Because the target genome is a single organelle contig (~15 kb animal mt,
~150 kb plant cp, up to ~1 Mb plant mt), the index is trivial to build
on the fly. We can omit `-d`/`-I` entirely and let miniprot index at
runtime — one-shot per sample.

## 6. What miniprot output tells us about assembly quality

Miniprot's per-hit features double as QC signals for C5:

| Feature in output | What it tells us | Action in C5 |
|---|---|---|
| Chain of exons with normal splice signals | Clean gene, well-assembled | Pass |
| Frameshift markers (`$` in alignment, `frameshift=` in GFF attributes) | Assembly error mid-CDS or genuine pseudogene | Fail the locus; flag sample in report |
| Internal stop codon (`*` in the `Prot:` line) | Wrong translation table *or* real premature stop | Try alternate table (animal mt), else fail |
| Truncated alignment (`Identity/Positive` low, short target span) | Divergent / partial hit | Length threshold check; fail if below `min_locus_length` |
| No hit at all | Locus missing (real gene loss or missed assembly) | Record in `validation.tsv` with `reason=no_hit` |
| Multiple non-overlapping hits (with `--outn > 1`) | Paralog / duplication | Log; retain primary only |

This is why we ask for `--gff --aln`: the alignment lines are what make
frameshift/stop-codon calls machine-readable rather than requiring us to
re-align externally.

## 7. What miniprot does *not* do (relevant to us)

- **Not a gene predictor.** It only reports alignments where an input
  protein hits — it cannot discover unannotated genes with no protein
  reference. This is fine for our purposes (we know the barcode panel in
  advance); it does mean the panel must be complete for a kingdom.
- **Not a splice-model expert.** The default splice model is a small
  general PWM (§2.9), less sophisticated than Spaln2's species-specific
  models. For animal mt (no introns) and plant plastid (few introns) this
  is a non-issue. For plant mt (introns present) it may cost a small
  amount of junction accuracy — acceptable for barcode extraction where
  we're aligning to the CDS, not annotating every exon-intron boundary.
- **Does not consume nucleotide references.** rRNA (12S, 16S) and tRNA
  (trnL) barcodes cannot be extracted with miniprot; those need a
  BLAST/nhmmer-based path if we later expand the panel beyond
  protein-coding loci ([tasks/3_refdata.md §2.3](../tasks/completed/3_refdata.md)
  notes this exclusion).
- **No profile/HMM support yet.** The paper flags HMMER-profile queries
  as future work — right now we hand it explicit protein FASTAs, which
  matches our `refs/.../proteins/<origin>/<gene>.faa` layout.

## 8. Integration notes for `MINIPROT_EXTRACT` (stage 13)

1. **Container:** biocontainer `quay.io/biocontainers/miniprot:0.18--h577a1d6_0`
   (pinned tag matches [conf/containers.config](../conf/containers.config)).
   ~5 MB image; no runtime dependencies beyond the binary.
2. **Input channel shape:** `tuple(meta, target_fasta) + path(locus_panel)` —
   `locus_panel` is the whole `proteins/` bundle root; the process picks
   the appropriate `<origin>/` subdirectory based on `meta.kingdom` +
   organelle inferred by `BIN_TARGET`.
3. **Runtime cap:** trivial — sub-minute per sample even with generous
   thread counts. Label as `process_low`.
4. **Determinism:** miniprot is deterministic given the same inputs and
   thread count; record miniprot version + panel bundle version in
   `metadata.json` for provenance.
5. **Failure mode:** if miniprot returns zero hits for a locus, treat as
   `no_recovery` for that locus (soft-fail in `validation.tsv`), **not**
   a Nextflow error. The pipeline still emits `barcodes.fasta` with
   whatever loci did recover, and the report surfaces the misses.
