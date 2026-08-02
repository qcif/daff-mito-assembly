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
| `METAFLYE` | `--genome-size` hint per target (§3.4). Flye's built-in polish (`--iterations 1`) runs identically for all targets. |
| `BANDAGE_NG` | No per-target variation (renders whatever graph is emitted). |
| `BIN_TARGET` | Per-target binning thresholds and criteria per §3.7. Homology is measured against the declared `${assembly_target}.mmi` **and** the sibling organelle panel(s), because plant samples carry both organelles and RECRUIT's positive selection does not separate them (§3.7.1). `plant_pt` runs the quadripartite canonicalisation of §3.6; `animal_mt` records circularity per §3.7.4. |
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
  IR resolution. For ORGANELLE_MAP visualisation it does; a concrete
  canonicalisation algorithm is specified in §3.6.

- **Plant mitogenome multi-isoform.** Plant mt recombines into multiple
  alternative arrangements; metaFlye emits multiple contigs with shared
  regions. Acceptable for barcoding — gene-bearing contigs are what matter.
  Diagnostics should flag when multiple alternative mt contigs are present.

- **Plastid carry-over into `plant_mt` assemblies.** A `plant_mt` run
  recruits against `plant_mt.mmi`, but plant plastid DNA is present at
  far higher copy number than mt DNA in the same extract and shares
  enough homology (plus genuine NUPT insertions in the mitogenome) to
  be recruited anyway. Observed on `INT-PLANT-01-mt` (§3.7): the two
  **highest-coverage** contigs in the assembly were plastid, aligning
  across 100 % of their length to `plant_pt.mmi` and < 7 % to
  `plant_mt.mmi`. Any binning rule that ranks on coverage alone will
  therefore pick the plastid as the mitogenome. Discrimination against
  the sibling organelle panel is mandatory on plant targets — see
  §3.7.1. The converse (mt carry-over into a `plant_pt` run) is far
  weaker, because the abundance gradient runs the other way.

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
  ~16 kb circle. Circularity is recorded by `BIN_TARGET` and is
  diagnostically valuable. **Flye's `circ.` column is the primary
  source, not end-overlap self-alignment** — see §3.7.4 for why the
  end-overlap check alone is a structural false negative on
  Flye-circularised contigs.

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


### 3.6 Plastid quadripartite canonicalisation

Plastid genomes have a canonical quadripartite structure: **LSC** (large
single-copy, ~80–90 kb) – **IRa** (inverted repeat A, ~20–30 kb) – **SSC**
(small single-copy, ~15–20 kb) – **IRb** (reverse-complement of IRa). Because
the IR is present twice in the plastome, a Flye assembly graph typically
emits **three edges**: LSC, SSC, and a single IR contig (whose read coverage
is roughly double the SC edges).

The two SSC orientations are biologically real — plastids exist as a mixture
of both isoforms — so a "correct" plastid assembly has two valid linear
representations (`path1` and `path2`). Our `BIN_TARGET` `plant_pt` branch
implements the algorithm described below as original in-house code;
see [spec/plastid-canonicalisation.md](plastid-canonicalisation.md) for
the full implementation specification.

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
5. On the canonical 3-edge branch, `BIN_TARGET` sets its primary `target.fasta` to `path1` and emits both isoforms as `plastid_isoforms/{path1,path2}.fasta`. **Precondition:** the substitution is conditional on C3 having selected at least one `plant_pt` contig (§3.7). Substituting unconditionally lets a sample whose assembly contains no recognisable plastid still emit a confident ~150 kb `target.fasta` — a silent false positive that defeats [principle 7](../CONSTITUTION.md) negative clarity. Where C3 selected nothing, C4 emits its isoforms as diagnostics only, `target.fasta` stays empty, and the disagreement is recorded in `bin_metadata.json` and surfaced in the report. Downstream stages (BLAST_VALIDATE, ANNOTATE, MINIPROT_EXTRACT) consume `target.fasta` unchanged — they are unaware of the plastid quadripartite structure and no plant_pt-specific branch exists inside them. `ORGANELLE_MAP` is the sole exception: when `plastid_isoforms/` is present it renders both isoforms. `bin_metadata.json` records the chosen LSC/IR/SSC edges, the canonicalisation branch, and both path lengths as the alternative-isoform provenance.

**Edge-count QC signal** is reported explicitly in the per-sample
`report.html` under the "Assembly quality assessment" section (§6a.2):
"3 edges → canonical", "1 edge → resolved circle", or "N edges → manual
review recommended". This gives the operator an immediate structural sanity
check without needing to open BandageNG.

Implementation: a clean-room stage-internal script
(`bin/plastid_canonicalise.py`) written from the algorithm description
above. Upstream is unmaintained and unlicensed, so re-implementing rather
than vendoring is the only viable option. The full algorithm
specification (input/output schema, decision tables, test matrix,
implementation-choice boundaries) lives in
[spec/plastid-canonicalisation.md](plastid-canonicalisation.md) — produced
by task 19 and consumed by task 20 as the sole permitted algorithm
reference during implementation.


### 3.7 Target binning criteria (per assembly target)

`BIN_TARGET` (C3) decides which assembled contigs constitute the declared
organelle. The original criterion ([task 18](../tasks/completed/18_bin_target.md))
was a flat, target-agnostic intersection:

> coverage spike (≥ 2× median contig coverage) **∩** reference identity
> (≥ 80 % **and** aligned fraction ≥ 0.30) **∩** ORF integrity (longest
> stop-free stretch ≥ 100 aa)

The first real `-profile integration` run to reach the assertion stage
(2026-07-31, all three fixtures) falsified three of its four assumptions.
Measured per-contig signals, with aligned fraction computed both as
implemented (single best alignment block) and as merged query intervals
across all blocks:

| Sample | Contig | bp | cov | best-block frac | **merged frac** | longest ORF (aa) | Truth |
|---|---|---|---|---|---|---|---|
| `INT-ANIMAL-01` | contig_8 | 16 952 | 71× | 0.477 | **0.878** | 201 | target (circular) |
| `INT-ANIMAL-01` | contig_9 | 16 942 | 31× | 0.733 | **0.985** | 337 | second mt copy |
| `INT-ANIMAL-01` | 13 others | 8.6–41 kb | 3–8× | ≤ 0.071 | **≤ 0.069** | 125–433 | off-target |
| `INT-PLANT-01-pt` | contig_2 | 136 976 | 148× | — | **0.998** | 606 | target |
| `INT-PLANT-01-pt` | contig_3 | 18 296 | 144× | — | **0.998** | 263 | target (SSC/IR) |
| `INT-PLANT-01-mt` | contig_1 | 69 287 | 168× | 0.016 | **0.064** | 570 | **plastid** (1.000 vs `plant_pt.mmi`) |
| `INT-PLANT-01-mt` | contig_6 | 86 004 | 137× | 0.007 | **0.012** | 773 | **plastid** (1.000 vs `plant_pt.mmi`) |
| `INT-PLANT-01-mt` | contig_3 | 30 195 | 7× | 0.128 | **0.353** | 254 | mt |
| `INT-PLANT-01-mt` | contig_5 | 84 652 | 5× | 0.082 | **0.165** | 369 | mt |

Net outcome under the original criterion: `animal_mt` selected correctly,
`plant_pt` and `plant_mt` selected **nothing at all**. The four findings
below are the spec response.

#### 3.7.1 Homology is measured against sibling panels too

`INT-PLANT-01-mt`'s two highest-coverage contigs are chloroplast (§3.3).
Recruitment cannot prevent this — it is positive selection against one
panel ([principle 5](../CONSTITUTION.md)), and plastid DNA is both far
more abundant and partly homologous (NUPTs). Separation is therefore
`BIN_TARGET`'s job, exactly as principle 5 intends.

**Change:** C3 loads the declared `${assembly_target}.mmi` *and* the
sibling organelle panel(s) from the same reference bundle, and a contig
is a candidate only where its merged aligned fraction against the
declared panel exceeds that against every sibling panel. Sibling sets:
`plant_mt` → `plant_pt`; `plant_pt` → `plant_mt`; `animal_mt` → none
required (no sibling organelle in the sample) but the check is run
uniformly for auditability.

This supersedes the previous §3.2 statement that C3 does "no cp/mt
discrimination inside the sample". That assumption held only while the
declared target was assumed to be the only organelle recruited; the
evidence above shows it is not.

#### 3.7.2 Aligned fraction is merged across all alignment blocks

The implemented metric took a single alignment block —
`best.blen / len(contig)`, where `best` is the block with the most
matching bases. minimap2 returns many local blocks per contig, so this
understates contig coverage by 2–10× on divergent input (0.477 → 0.878
on the animal target; 0.128 → 0.353 on the plant mt target). This is a
defect in the metric, not a threshold: it does not measure what its name
says.

**Change:** aligned fraction is the union of query intervals over all
alignment blocks against a panel, divided by contig length. Identity is
reported as `sum(matching bases) / sum(block lengths)` over that same
merged set rather than one block's identity.

On the corrected metric the signal is cleanly bimodal where it should be
— `animal_mt` separates 0.878/0.985 (target) from ≤ 0.069 (off-target),
and both `plant_pt` contigs sit at 0.998. **`plant_pt` needs no threshold
relaxation whatsoever**; it failed solely on this defect plus §3.7.3.

#### 3.7.3 Coverage is a ranking signal, not a gate

The coverage-spike gate assumed a nuclear background to spike above.
METAFLYE consumes RECRUIT output, so that background has already been
removed — the median contig coverage in an organelle assembly *is*
target coverage. The gate is consequently self-defeating: it passes only
when recruitment leaves enough junk behind to depress the median
(`animal_mt`, 13 low-coverage survivors → median 4×) and fails when
recruitment works well (`plant_pt`, 2 clean contigs → median 146×,
requiring an impossible 292×).

Worse, on `plant_mt` it is **inverted**: the only contig clearing
`2× median` was contig_1 — a plastid. A relaxation of the identity
threshold alone, with the coverage gate retained, would have emitted
155 kb of chloroplast as the plant mitogenome.

**Change:** coverage ceases to be an admission gate. It remains the
**dominance** signal used to rank admitted candidates (`animal_mt` and
`plant_pt` emit the single highest-coverage candidate) and is retained
in `secondaries.tsv` and the report as a diagnostic.

**Retained risk — NUMTs/NUPTs.** The coverage gate was the incidental
defence against binning a nuclear insertion of organellar DNA, which is
homologous but at nuclear copy number. Under the revised criteria that
defence comes from §3.7.1 plus dominance ranking, not from an absolute
coverage floor. This is a deliberate trade: a per-target absolute floor
cannot be set without excluding genuine low-coverage plant mt
sub-genomic contigs (5–7× in the fixture above, against a 137–168×
plastid in the same assembly). NUMT risk is therefore accepted and
**flagged** rather than filtered: a candidate whose coverage falls below
`low_coverage_fraction` (0.05) of the top-ranked candidate's is still
emitted, but carries `low_coverage_candidate: true` in
`secondaries.tsv` and `bin_metadata.json` for the report to surface
([principle 7](../CONSTITUTION.md) — the operator sees the uncertainty
rather than the pipeline silently guessing). See
[task 23](../tasks/completed/23_bin_target_recalibration.md) §5.2 for the
flag, and [task 26](../tasks/26_binning_marker_genes.md) for the
marker-gene criterion that would separate NUMTs properly if the flag
proves insufficient.

#### 3.7.4 Circularity comes from Flye, not from end-overlap

Flye **trims the terminal overlap** when it circularises a contig, so a
Flye-circularised contig has no residual end redundancy left to detect.
Self-aligning the first and last *N* bp of `INT-ANIMAL-01`'s target
returns no hit at N = 300, 1 000, or 5 000 — while `assembly_info.txt`
records `circ. = Y` for that contig. The end-overlap check is a false
negative by construction on exactly the contigs it exists to confirm; it
can only ever fire on a contig Flye *failed* to circularise.

**Change:** the `circ.` column of `assembly_info.txt` is the primary
signal. The end-overlap self-alignment is retained as a fallback for
contigs Flye marks `N`. `bin_metadata.json` records which method fired
(`flye_circ` | `end_overlap` | `none`), because Flye's own circularity
call can be over-confident on repeat-collapsed contigs and the auditor
needs to see which evidence was used ([rule 16](../CONSTITUTION.md)).

#### 3.7.5 ORF integrity as implemented is vacuous

The longest stop-free stretch exceeds the 100 aa floor for **every
contig in every fixture** (125–773 aa), including an 8.6 kb contig with
zero reference homology. At these contig lengths a 100 aa stop-free
stretch arises by chance in any sequence; the criterion contributes no
discrimination, and in the `plant_mt` fixture it is anti-correlated with
truth (the longest ORF, 773 aa, sits on a plastid contig).

**Change:** the longest-ORF measure is demoted to a recorded diagnostic
and drops out of the selection intersection. It is *not* replaced by a
raised threshold — the measure is the wrong one, not badly calibrated.
The principled replacement is **panel marker-gene presence** (does the
contig carry target-appropriate organellar protein-coding genes?), which
would also give `plant_mt` a positive discriminator. That change needs a
gene set broader than the `plant_mt` barcode panel of §4.3, and
therefore a reference-bundle addition; it is specified as a conditional
change in [task 26](../tasks/26_binning_marker_genes.md), gated on
whether the §3.7.1 sibling test proves sufficient on its own, rather
than adopted here.

#### 3.7.6 Revised per-target criteria

A contig is a **target candidate** where all of:

1. merged aligned fraction against the declared panel ≥ `min_aligned_frac`;
2. merged identity against the declared panel ≥ `min_identity`;
3. merged aligned fraction against the declared panel > that against
   every sibling panel (§3.7.1).

| `assembly_target` | `min_identity` | `min_aligned_frac` | Sibling panel | Emit | Rationale |
|---|---|---|---|---|---|
| `animal_mt` | 80 % | 0.30 | — | single highest-coverage candidate | Observed separation 0.878/0.985 vs ≤ 0.069 — wide margin either side of 0.30. |
| `plant_pt`  | 75 % | 0.30 | `plant_mt` | single highest-coverage candidate | Both true contigs at 0.998; the 75 % floor accommodates the SSC/IR contig at 79.1 %. C4 substitution then applies per §3.6. |
| `plant_mt`  | 70 % | 0.15 | `plant_pt` | all candidates, up to 20 contigs | Plant mt is ~90 % fast-evolving non-coding sequence, so whole-genome homology to a non-conspecific reference is genuinely low (0.165–0.353 for true mt). The sibling-panel test in §3.7.1, not the floor, is what excludes the plastid. |

These are **prototype defaults calibrated on three fixtures**, not
validated thresholds. They remain subject to
[§9 item 10](07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
benchmarking, and per [rule 18](../CONSTITUTION.md) they live in
versioned config rather than as constants in `bin/bin_target.py` — as
`params.bin_target_thresholds` in `nextflow.config`, keyed by
`assembly_target`, passed to C3 as CLI arguments and echoed back into
`bin_metadata.json` as `thresholds_applied`. Each entry also carries
`emit` (`single` | `all`), `max_contigs`, and the
`low_coverage_fraction` of §3.7.3.

Measured on the 2026-08-02 `-profile integration` run under these
values (task 23), no threshold required adjustment: `INT-ANIMAL-01`
selects contig_8 (circular via `flye_circ`), `INT-PLANT-01-pt` selects
contig_2 on C3's own evidence, and `INT-PLANT-01-mt` selects contig_3 +
contig_5 (114 847 bp) while classifying contig_1 and contig_6 as
`sibling_organelle` at 1.000 against `plant_pt.mmi`.

**`plant_mt` residual uncertainty.** Under these criteria the
`INT-PLANT-01-mt` fixture yields contig_3 + contig_5 ≈ 115 kb, below the
200 kb – several Mb range quoted in §3.1. Either the mitogenome is
genuinely incompletely assembled at 5–7× coverage, or contigs are being
missed. The former `expected/plant_mt/bin_bounds.json` bound of
200 000–700 000 bp was **not** evidence against the criteria — it was
derived from the whole-assembly total, which is now known to be 57 %
plastid, and treating it as a target to hit would mean binning
chloroplast as mitochondrion. Task 23 corrected it to a deliberately
wide 90 000–400 000 bp. It should not be tightened until
[task 25](../tasks/25_coverage_gate_carryover.md) establishes whether
the recruited read pool, not C3, is what limits the assembly.
