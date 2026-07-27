## 3. Organelle considerations: mitochondrion vs plastid

Plants carry both a mitogenome and a plastome; animals carry only a
mitogenome. The two are similar enough that a single assembly pass handles
both, but several downstream parameters must vary by organelle type. This
section enumerates what differs and where the workflow forks.

### 3.1 What differs between organelles

| Aspect | plant cp | plant mt | animal mt |
|---|---|---|---|
| Typical size | 120–160 kb | 200 kb – several Mb | 15–20 kb |
| Structure | quadripartite, large inverted repeat | multipartite, recombinant, multiple isoforms | simple circle, compact |
| Relative coverage | very high (>1000× possible) | medium | medium–high |
| Genetic code (ORF validation) | NCBI table 11 (bacterial/plastid) | NCBI table 1 (standard) | table 2 (vertebrate) / table 5 (invertebrate) — clade-dependent |
| Primary barcodes | rbcL, matK, trnH-psbA, ndhF, atpB | rarely used for barcoding | COX1, CYTB, 12S, 16S |
| Annotation tool | TBD ([§8](07-open-questions.md#8-remaining-open-questions)) | TBD ([§8](07-open-questions.md#8-remaining-open-questions)) | MITOS2 (or MitoFinder) |
| BLAST validation DB | RefSeq plastid | RefSeq plant mitochondrion | RefSeq metazoan mitochondrion |
| Assembly pitfall | inverted repeat collapses or branches ambiguously — canonical quadripartite form not guaranteed (canonicalised post-assembly per §3.6) | metaFlye emits multiple alternative isoforms; no canonical form | should close cleanly as one circle; verify circularity |

### 3.2 Fork vs dynamic parameter — by stage

| Stage | Approach | Rationale |
|---|---|---|
| `RECRUIT` | **dynamic**: kingdom-keyed reference panel bundles all relevant organelles (plant = cp + mt in one panel) | A single recruitment pass keeps the read pool unified for assembly. Plant kingdom panel must include both cp and mt references. |
| `METAFLYE` | **shared**: one assembly pass per sample | No assembly-time fork. metaFlye handles mixed-organelle input fine; contigs separate naturally by coverage and gene content. |
| `MEDAKA` | **shared, opt-in**: one polish pass when `--polish` is passed | Identical model regardless of organelle. Skipped by default; enable with `--polish`. |
| `BANDAGE_NG` | **shared**: one graph image | The full graph shows cp + mt branching for plants — diagnostically useful. |
| `BIN_TARGET` | **fork point**: classify each contig as cp / mt / off-target by reference identity + coverage tier | Natural fork point. Downstream stages branch per-contig from here. |
| `BLAST_VALIDATE` | **dynamic**: per-organelle BLAST DB selected by `BIN_TARGET` label | Validation must match the contig's classified organelle. |
| `MINIPROT_EXTRACT` | **dynamic**: per-contig genetic code table and locus panel subset selected by organelle label (and clade-trial for animal mt — see §3.3) | Genetic code differs by organelle; locus panel partitions naturally (rbcL → cp, COX1 → mt). |
| `ANNOTATE` | **dynamic**: MITOS2 for animal mt, TBD for plant cp+mt ([§8](07-open-questions.md#8-remaining-open-questions)) (per organelle label from `BIN_TARGET`) | Annotator is organelle-specific. |
| `ORGANELLE_MAP` | **shared**: consumes annotation output from `ANNOTATE` | Renderer is generic over organelle type. |
| `COLLATE` | **shared**: aggregate all per-organelle outputs into one bundle | Plant samples produce cp + mt barcode FASTAs concatenated into a single emission. |

### 3.3 Specific issues and decisions

- **Plant samples emit two organelles per target organism.** Animal samples
  emit one (mt); plant samples emit both cp and mt of the same target plant.
  Brief.md §3.5 wording covers this via "dominant kingdom-matched organelle
  assembly" — for plants both organelles of the dominant organism are in
  scope.

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

- **Animal mt genetic code is clade-dependent.** Vertebrate (table 2),
  invertebrate (table 5), echinoderm (table 9), etc. The submitter declares
  kingdom but not phylum. Options: (a) default to invertebrate (table 5) —
  covers most biosecurity intercepts; (b) require phylum at intake; (c) try
  table 2 and table 5 in EXTRACT, pick the one with valid ORF.
  **Recommend (c)** — automatic, unambiguous when an ORF is recoverable, and
  avoids burdening the submitter.

- **rDNA out of scope** ([brief.md §2](brief.md)). Nuclear rDNA (ITS, 18S,
  28S) recovery is dropped from this workflow. Extraction is limited to
  organelle-encoded loci.

- **Animal mt circularisation check.** Animal mt should close as a single
  ~16 kb circle. End-overlap circularity check added to `BIN_TARGET` is
  diagnostically valuable; flag in P3.

### 3.4 Flye `--genome-size` hint (per organelle)

Flye accepts a `--genome-size` hint (e.g. `--genome-size 135k`) that improves
assembly on small, high-coverage targets like organelles.
[CLAW](../CLAW/config.yml) hard-codes 135 kb for chloroplasts; we generalise per
organelle label emitted by `BIN_TARGET` — except Flye runs *before*
`BIN_TARGET`, so the hint at the METAFLYE stage is set per-**kingdom** using
the largest expected organelle in that kingdom (so we don't under-hint a
plant mitogenome).

| Kingdom | Hint at METAFLYE (largest expected) | Rationale |
|---|---|---|
| plant | `2m` | Plant mt can reach 200 kb–several Mb; cp is small enough not to be prejudiced. |
| animal | `20k` | Animal mt is ~15–20 kb; tight hint helps recover the small circle. |

The hint is advisory to Flye's coverage estimator, not a length filter — over-
or under-hinting degrades but does not silently truncate. Re-evaluate the
plant hint in P2 if mt assembly fragments.

### 3.5 Flye `--asm-coverage` for plastid runs

For plant samples, when a plastid contig is being targeted at high coverage,
pass Flye `--asm-coverage <N>` where `N = coverage_max − 20` (COVERAGE_GATE's
plant MAX minus 20; see §2.1.1). Rationale (adopted from
[ptGAUL](../reference-material/ptgaul/ptGAUL.sh)): when total coverage exceeds
this value, Flye uses only the longest reads for the initial contig-building
step, which improves contiguity on high-copy small circles like plastids.
Skip for animal mt — the target is too small for `--asm-coverage` to matter.

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
lifted directly into our `BIN_TARGET` plant-cp branch.

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
