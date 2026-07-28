## 3. Organelle considerations: mitochondrion vs plastid

Plants carry both a mitogenome and a plastome; animals carry only a
mitogenome. Per [§1a](01-pipeline-flow.md#1a-engineering-constraints),
each sample-sheet row targets exactly one organelle — `plant_pt`,
`plant_mt`, or `animal_mt` — so an operator wanting both plastid and
mitogenome from a single plant sample submits two rows sharing the same
reads. This keeps every stage of the pipeline monomorphic per run:
one reference index, one coverage limit, one genome-size hint, one
genetic-code table. The differences enumerated below therefore select
by `meta.assembly_target` (a scalar) rather than fanning out
per-organelle within a single sample.

### 3.1 What differs between assembly targets

| Aspect | `plant_pt` | `plant_mt` | `animal_mt` |
|---|---|---|---|
| Typical size | 120–160 kb | 200 kb – several Mb | 15–20 kb |
| Structure | quadripartite, large inverted repeat | multipartite, recombinant, multiple isoforms | simple circle, compact |
| Relative coverage | very high (>1000× possible) | medium | medium–high |
| Genetic code (ORF validation) | NCBI table 11 (bacterial/plastid) | NCBI table 1 (standard) | table 2 (vertebrate) / table 5 (invertebrate) — clade-dependent |
| Primary barcodes | rbcL, matK, trnH-psbA, ndhF, atpB | rarely used for barcoding | COX1, CYTB, 12S, 16S |
| Annotation tool | TBD ([§8](07-open-questions.md#8-remaining-open-questions)) | TBD ([§8](07-open-questions.md#8-remaining-open-questions)) | MITOS2 (or MitoFinder) |
| BLAST validation DB | RefSeq plastid | RefSeq plant mitochondrion | RefSeq metazoan mitochondrion |
| Assembly pitfall | inverted repeat collapses or branches ambiguously — canonical quadripartite form not guaranteed (canonicalised post-assembly per §3.6) | metaFlye emits multiple alternative isoforms; no canonical form | should close cleanly as one circle; verify circularity |

### 3.2 Per-stage parameter selection

Every stage keys off `meta.assembly_target` (scalar per-sample); no
stage forks *within* a single sample. The tables in §3.1 and the
per-target coverage limits ([spec §2.1.1](02-stages.md#211-per-target-limits))
supply the concrete values.

| Stage | Target-driven parameter(s) |
|---|---|
| `RECRUIT` | Single reference index: `${assembly_target}.mmi` (one minimap2 pass). |
| `METAFLYE` | `--genome-size` hint per target (§3.4); `--asm-coverage` on `plant_pt` only (§3.5). |
| `MEDAKA` | No per-target variation (identical model). Opt-in via `--polish`. |
| `BANDAGE_NG` | No per-target variation (renders whatever graph is emitted). |
| `BIN_TARGET` | Reference identity check against the same `${assembly_target}.mmi` used at RECRUIT; only classifies contigs as target vs off-target — no cp/mt discrimination inside the sample. `plant_pt` runs the quadripartite canonicalisation of §3.6; `animal_mt` runs the end-overlap circularity check (§3.3). |
| `BLAST_VALIDATE` | BLAST DB fixed per-target: `refseq_pt` for `plant_pt`, `refseq_mt_viridiplantae` for `plant_mt`, `refseq_mt_metazoa` for `animal_mt`. |
| `MINIPROT_EXTRACT` | Locus panel subset per target (rbcL/matK for `plant_pt`; cox1/nad1 for `plant_mt`; COX1/CYTB for `animal_mt`). Genetic-code table per target (11 / 1 / 2-or-5). Clade-trial only on `animal_mt` (§3.3). |
| `ANNOTATE` | Annotator per target — MITOS2 for `animal_mt`; TBD ([§8](07-open-questions.md#8-remaining-open-questions)) for `plant_pt` / `plant_mt`. |
| `ORGANELLE_MAP` | No per-target variation (renderer is generic over annotated GenBank input). |
| `COLLATE` | One per-sample bundle per row — no cross-target aggregation inside a sample. |

### 3.3 Specific issues and decisions

- **One sample row → one organelle bundle.** Animal samples submit a
  single `animal_mt` row and get one bundle. Plant samples requiring
  both organelles submit two rows (`plant_pt` + `plant_mt`) with the
  same reads; each produces its own independent bundle under a distinct
  `sample_id`. Brief.md §3.5 wording ("dominant kingdom-matched organelle
  assembly") is now per-target: the dominant contig **for the declared
  `assembly_target`** is emitted; secondaries are recorded as
  diagnostics.

- **Plastid inverted repeat.** metaFlye may collapse the IR or emit alternative
  paths; canonical quadripartite form (LSC–IRb–SSC–IRa) is not guaranteed. For
  barcode extraction this does not matter — miniprot finds genes regardless of
  IR resolution. For ORGANELLE_MAP visualisation it does; a concrete canonicalisation
  algorithm derived from [ptGAUL](../reference-material/ptgaul/ptGAUL.sh) is
  specified in §3.6.

- **Plant mitogenome multi-isoform.** Plant mt recombines into multiple
  alternative arrangements; metaFlye emits multiple contigs with shared
  regions. Acceptable for barcoding — gene-bearing contigs are what matter.
  Diagnostics should flag when multiple alternative mt contigs are present.

- **`animal_mt` genetic code is clade-dependent.** Vertebrate (table 2),
  invertebrate (table 5), echinoderm (table 9), etc. The submitter declares
  the assembly target but not the phylum. Options: (a) default to
  invertebrate (table 5) — covers most biosecurity intercepts; (b) require
  phylum at intake; (c) try table 2 and table 5 in EXTRACT, pick the one
  with valid ORF. **Recommend (c)** — automatic, unambiguous when an ORF
  is recoverable, and avoids burdening the submitter.

- **rDNA out of scope** ([brief.md §2](brief.md)). Nuclear rDNA (ITS, 18S,
  28S) recovery is dropped from this workflow. Extraction is limited to
  organelle-encoded loci.

- **`animal_mt` circularisation check.** Animal mt should close as a single
  ~16 kb circle. End-overlap circularity check added to `BIN_TARGET` is
  diagnostically valuable; flag in P3.

### 3.4 Flye `--genome-size` hint (per target)

Flye accepts a `--genome-size` hint (e.g. `--genome-size 135k`) that improves
assembly on small, high-coverage targets like organelles.
[CLAW](../CLAW/config.yml) hard-codes 135 kb for chloroplasts; we set the
hint per-target directly from `meta.assembly_target`, since each run
assembles one organelle.

| Assembly target | METAFLYE `--genome-size` | Rationale |
|---|---|---|
| `plant_pt` | `150k` | Land-plant plastids sit tightly around 120–160 kb. |
| `plant_mt` | `2m` | Plant mt can reach 200 kb–several Mb; upper bound keeps Flye's coverage estimator honest. |
| `animal_mt` | `20k` | Animal mt is ~15–20 kb; tight hint helps recover the small circle. |

The hint is advisory to Flye's coverage estimator, not a length filter —
over- or under-hinting degrades but does not silently truncate.
Re-evaluate the `plant_mt` hint in P2 if mt assembly fragments.

### 3.5 Flye `--asm-coverage` for `plant_pt` runs

For `plant_pt` samples, pass Flye `--asm-coverage <N>` where
`N = coverage_max − 20` (`plant_pt` MAX from §2.1.1 minus 20).
Rationale (adopted from [ptGAUL](../reference-material/ptgaul/ptGAUL.sh)):
when total coverage exceeds this value, Flye uses only the longest reads
for the initial contig-building step, which improves contiguity on
high-copy small circles like plastids. Skipped for `plant_mt` and
`animal_mt` — the plant mt is too large to be constrained this way, and
the animal mt too small for `--asm-coverage` to matter.

### 3.6 Plastid quadripartite canonicalisation (ptGAUL-derived)

Plastid genomes have a canonical quadripartite structure: **LSC** (large
single-copy, ~80–90 kb) – **IRa** (inverted repeat A, ~20–30 kb) – **SSC**
(small single-copy, ~15–20 kb) – **IRb** (reverse-complement of IRa). Because
the IR is present twice in the plastome, a Flye assembly graph typically
emits **three edges**: LSC, SSC, and a single IR contig (whose read coverage
is roughly double the SC edges).

The two SSC orientations are biologically real — plastids exist as a mixture
of both isoforms — so a "correct" plastid assembly has two valid linear
representations (`path1` and `path2`). This is implemented in
[ptGAUL's combine_gfa.py](../reference-material/ptgaul/combine_gfa.py) and is
lifted directly into our `BIN_TARGET` `plant_pt` branch.

**Algorithm:**

1. Parse Flye's `assembly_graph.gfa`; extract per-edge sequence (from `S` lines) and per-edge depth (from the `S` line depth tag).
2. Count edges. Classify:
   - **3 edges** → canonical quadripartite. Proceed to canonicalise.
   - **1 edge** → fully-resolved circle. Passthrough as the assembly; flag "no IR resolution needed" in the report.
   - **anything else** → non-canonical. Emit as-is with a diagnostic flag; do not attempt automated canonicalisation. Human review via BANDAGE_NG output.
3. For the 3-edge case:
   - Longest edge by sequence length → **LSC**.
   - Deepest edge by read depth (top of `sort -k2 -n -r`) → **IR** (2× coverage betrays it).
   - Remaining edge → **SSC**.
4. Emit two canonical isoforms:
   - `path1.fasta`: `LSC + IR + SSC + reverse_complement(IR)`
   - `path2.fasta`: `LSC + IR + reverse_complement(SSC) + reverse_complement(IR)`
5. Record both isoforms in per-sample outputs. Downstream stages (BLAST_VALIDATE, ANNOTATE, ORGANELLE_MAP, MINIPROT_EXTRACT) run against `path1` as the primary; `path2` is emitted alongside and noted in `metadata.json` as the alternative isoform.

**Edge-count QC signal** is reported explicitly in the per-sample
`report.html` under the "Assembly quality assessment" section (§6a.2):
"3 edges → canonical", "1 edge → resolved circle", or "N edges → manual
review recommended". This gives the operator an immediate structural sanity
check without needing to open BandageNG.

Implementation: port `combine_gfa.py` (small, ~90-line Python) into the
pipeline as a stage-internal script; upstream is unmaintained and pinning it
directly is cleaner than a conda dep.
