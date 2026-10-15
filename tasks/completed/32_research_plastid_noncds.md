# Task 32 — Research: plastid (`plant_pt`) non-CDS annotation

**Phase:** P4 (from [spec §6](../spec/06-phases.md)), **research only —
no implementation.** Produces a decision and its evidence, not code.

**Status:** done (2026-10-15). Decision implemented by task 65_plastid_noncds_barrnap.md. The Chloë licence enquiry (§7.3) was dropped by the user, so Chloë is not pursued.

## 0. Overview

[Task 30](completed/30_unified_locus_pass.md) gives every target a protein-coding
annotation from the unified miniprot pass, and
[task 31](completed/31_annotate_merge.md) publishes it. A `plant_pt` sample
therefore already receives a real annotation of its ~79 plastid
protein-coding genes, with `status: "ok_cds_only"`.

What it does **not** receive:

- **tRNA features** (~30 in a typical plastome);
- **rRNA features** (4 distinct, 8 loci counting the inverted repeat);
- **correct intron and *trans*-splicing structure** — `rps12` is split
  across the genome and reassembled in *trans*, and several plastid
  genes carry group II introns. miniprot reports spliced alignments but
  will call these as fragments or miss them, so the CDS model for those
  specific genes is incomplete.

This task decides whether to fill those gaps, and with what.

**Why this is no longer urgent.** Before the unified locus pass, this
question blocked P4 — a plant sample would have received an *empty*
annotation. It now receives a good one that is explicitly labelled
CDS-only. The remaining question is an enhancement, weighed against
[rule 19](../CONSTITUTION.md): is a complete plastid feature map worth a
new dependency, given that stage 13a is supplementary to the contractual
barcode path ([spec §2 stage 13](../spec/02-stages.md#2-stage-detail))?

**Nothing downstream is blocked.** `ok_cds_only` keeps the channel shape
uniform, so `COLLATE` and `RUN_REPORT` proceed normally. P4 completes
without this task.

## 1. Evaluated and ruled out (with reasons)

Recorded so the search is not repeated. Verified 2026-08-04.

| Candidate | Status | Reason |
|---|---|---|
| **GeSeq** | Ruled out | Web service, no CLI. A containerised pipeline cannot use it ([rule 11](../CONSTITUTION.md)). Named in [brief.md §3.6](../brief.md); the brief should be reconciled. |
| **Chloë** ([`ian-small/Chloe.jl`](https://github.com/ian-small/Chloe.jl)) | **Blocked, not rejected** | Best technical fit — purpose-built by the Small lab (UWA), GFF3 by default, handles introns and the IR, XGBoost-scored features. **But: no LICENSE file on either the code or the [references repo](https://github.com/ian-small/chloe_references)** (`LICENSE` 404s on both; GitHub API reports no licence), so default copyright applies and we may not redistribute it in an image. **And:** the README states it is optimised for angiosperms and directs non-angiosperm users to contact the author — a real limit for a workflow whose premise is that the submitter's taxon is unknown ([brief.md §2](../brief.md)); imported timber and conifer nursery stock are routine Australian interceptions. Also Julia (bespoke Dockerfile, [rule 12](../CONSTITUTION.md) last resort), no git tags, 5 stars, last push 2025-06-09. |
| **PGA** ([`quxiaojian/PGA`](https://github.com/quxiaojian/PGA)) | Weak fallback | GPL-3.0 — licensing clean. But Perl, last pushed 2020-10-29, requires user-supplied reference GenBanks. Licensed-but-stale. |
| **CPGAVAS2** | Not evaluated | Web + heavyweight Docker (MySQL-backed). Evaluate only if everything above closes out. |
| **MITOS2** | Ruled out | Metazoan mitochondrial models and data. On a plastome it would produce confident nonsense. |

**No bioconda recipe exists for any plastid annotator** — checked
directly against the quay.io biocontainers API. That is why
[rule 12](../CONSTITUTION.md)'s bespoke-Dockerfile last resort is on the
table at all.

## 2. The cheap option, evaluated first

Before reaching for a plastid-specific annotator, measure what
**tRNAscan-SE + barrnap** give. Both are bioconda
(`trnascan-se` 2.0.13, `barrnap` 1.10.6), both drop into
[task 31](completed/31_annotate_merge.md)'s existing non-CDS slot with no new
architecture, and together they address two of the three gaps:

- `tRNAscan-SE -O` (organellar mode) for tRNAs;
- `barrnap` for rRNAs — plastid rRNAs (`rrn16`, `rrn23`, `rrn5`,
  `rrn4.5`) are bacterial-like, so the **bacterial** models are the
  right ones. Verify empirically; barrnap's `mito` models are metazoan
  and are the wrong choice here.

They do **not** address intron structure. If the measurement shows
tRNA/rRNA recovery close to the reference annotation, the residual gap
is `rps12` and a handful of intron-containing genes — which may well be
acceptable for a diagnostic map, and should then be *recorded* in the
summary rather than fixed.

**This is the recommended first experiment**, because it is a
`mulled-build` of two bioconda packages
([spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)
option 2) rather than a bespoke Julia image with a licence problem.

## 3. Measurement

**Fixture and ground truth.** `INT-PLANT-01-pt`'s binned `target.fasta`
is the canonicalised *Datura stramonium* plastome, whose top
`BLAST_VALIDATE` hit is `NC_018117.1` at 99.60 % identity
([task 21 §9](completed/21_blast_validate.md#9-outcomes)). Use that
RefSeq record's own annotation as truth. This is a legitimate use of a
reference annotation — measuring an annotator against a curated answer,
not tuning a threshold against a fixture bound
([spec §5.1](../spec/05-test-data.md)).

Measure, for §2's option and for Chloë if §4 unblocks it:

1. tRNAs recovered vs the reference (~30) and spurious calls;
2. rRNAs recovered (4 distinct / 8 loci with the IR) — **and whether
   both IR copies are called**;
3. intron-containing genes and `rps12` handled correctly;
4. wall time and peak RSS against the 60-minute CI budget.

**Accept §2 if** tRNA recovery ≥ 90 % with < 10 % spurious, and both IR
rRNA copies are called. Then the only recorded limitation is intron
structure, and no bespoke container is needed.

## 4. Open questions

1. **Chloë licence.** Email Ian Small (`ian.small@uwa.edu.au`) and Ian
   Castleden — both named in `Project.toml`, and the README invites
   contact — requesting an explicit open-source grant (MIT/Apache-2.0)
   covering code and references, or written permission for DAFF to
   redistribute it in a container image. **Likely to succeed** —
   publicly-funded Australian academic software, free public web
   service, Australian government requester; the missing licence reads
   as oversight. Cheap to ask, so ask early even though §2 may make it
   moot. Precedent: task 19/20 hit the same problem with the upstream
   canonicalisation code and resolved it by re-implementing, which is
   not viable at Chloë's scale.
2. **Chloë non-angiosperm coverage.** Ask in the same email.
3. **Does anyone need it?** Establish from P3/P5 testing or from the
   client whether an operator needs a complete plastid feature map, or
   whether the CDS annotation plus a stated limitation suffices. This
   decides between §2, Chloë, and doing nothing — and cannot be answered
   from the code.

## 5. Exit criteria

- §2 measured against §3's criteria.
- Chloë licence and scope questions asked and answered (or a recorded
  non-response after a reasonable interval).
- A decision recorded in
  [spec §3.1/§3.2](../spec/03-organelles.md#31-what-differs-between-assembly-targets)
  with its evidence and rejected alternatives, so an auditor sees why
  the plastid annotation is what it is without opening a task file
  ([rule 18](../CONSTITUTION.md)).
- [spec §8 item 3](../spec/07-open-questions.md#8-remaining-open-questions)
  closed for `plant_pt`.
- An implementation task written, **or** a recorded decision that
  `ok_cds_only` stands — a legitimate outcome, not a failure to decide.

## 6. Notes

- **Do not fork or vendor Chloë** pending §4 item 1, and do not build an
  image on the assumption a grant will arrive.
- **`plant_mt` non-CDS is [task 37](shelved/archived/37_annotate_plant_mt_noncds.md)**,
  blocked on its own preconditions. If §2 is accepted here, task 37
  becomes nearly free — the same two tools serve both plant targets.

## 7. Outcomes

### 7.1 Measurement (2026-10-13) — §2 option

**Setup.** `INT-PLANT-01-pt` `bin_target/target.fasta` (one contig,
155,276 bp, canonicalised), against `NC_018117.1` (155,871 bp, fetched
from NCBI efetch). Reference tRNA/rRNA features were lifted onto the
fixture by `blastn` of each feature's full genomic span (introns
included, which is also how tRNAscan-SE reports a spliced tRNA),
keeping hits ≥ 90 % identity over ≥ 90 % of the query. A call matches
a reference locus if it covers ≥ 50 % of it and has the same anticodon
(tRNA) or the same rRNA species. Scratch scripts were not committed.

Containers (both bioconda, both checked with `--help`):

- `quay.io/biocontainers/trnascan-se:2.0.13--pl5321hab16a5f_0`, run
  as `tRNAscan-SE -O`;
- `quay.io/biocontainers/barrnap:1.10.6--pl5321hdfd78af_0`, run as
  `barrnap --kingdom bac`. Note: this image's `--version` reports
  **1.10.5**. barrnap 1.10 is a rewrite on Infernal + Rfam: it has no
  `mito`/`euk` kingdoms at all (only `bac`/`arc`/`fun`), and it can
  also call tRNAs (`--trna`, Rfam RF00005). `--all` **crashes** in
  this image (a broken `transterm` data path), so only `--rrna` and
  `--trna` are usable.

**Reference corrections.** `NC_018117.1` has curation errors that
would otherwise distort the score:

- `trnG-GCC` / `trnG-UCC` are swapped. The record puts `trnG-GCC` on
  the intron-containing locus (~9.3 kb) and `trnG-UCC` on the
  single-exon locus (~37.9 kb). Canonical angiosperm plastomes are
  the reverse, and both tools read **GCC** from the sequence at
  37.9 kb. Scored with the corrected labels.
- `trnP-GGG` (68752–68823) is a second annotation of the `trnP-UGG`
  locus (68750–68824), not a separate gene. Dropped.
- Several `product` fields are wrong (`trnH-GUG` "tRNA-Gly",
  `trnP-UGG` "tRNA-Trp", `trnI-CAU` "tRNA-Met"). Scored on `gene`.

That leaves **37 tRNA loci** and **8 rRNA loci** (4 species × 2 IR
copies).

**tRNA results** (both tools give the same full calls):

| Class | Loci | tRNAscan-SE `-O` | barrnap `--trna` |
|---|---|---|---|
| Intronless tRNAs | 29 | **29/29** | **29/29** |
| Intron-containing tRNAs (trnK-UUU, trnG-UCC, trnL-UAA, trnV-UAC, trnI-GAU ×2, trnA-UGC ×2) | 8 | 0 correct. 4 partial (trnI-GAU ×2, trnA-UGC ×2: right locus and anticodon, truncated structure, score ~16 against 56–79 for real calls). 4 missed. | 0 (all missed) |
| **Recovery (strict)** | 37 | **29/37 = 78.4 %** | **29/37 = 78.4 %** |
| Recovery (locus found, any structure) | 37 | 33/37 = 89.2 % | 29/37 = 78.4 % |
| Spurious calls | — | **0** | **0** |

Every call lies on a real locus. The three CAU tRNAs (`trnfM`, `trnM`,
`trnI-CAU` ×2) are all named Met by both tools. That is an isotype
ambiguity: the anticodon is correct but the amino acid is wrong for
`trnI-CAU`.

**rRNA results** (barrnap `--kingdom bac`):

| rRNA | IR copy A | IR copy B |
|---|---|---|
| rrn16 | ✓ | ✓ |
| rrn23 | ✓ (3′ end over-extended by ~94 bp) | ✓ (same) |
| rrn5 | ✓ | ✓ |
| rrn4.5 | ✗ | ✗ |

6/8 loci. **Both IR copies are called for every rRNA species that is
found.** rrn4.5 is missed in both copies; Rfam's bacterial LSU model
has no separate 4.5S, which is the plastid-specific split-off 3′ end
of 23S.

**Cost** (`/usr/bin/time -v` inside the container, 1 thread):

| Tool | Wall | Peak RSS |
|---|---|---|
| tRNAscan-SE `-O` | 86 s | 151 MB |
| barrnap `--rrna` | 6 s | **3.2 GB** |
| barrnap `--trna` | 6 s | **3.2 GB** |

Both fit inside the 60-minute CI budget. barrnap's 3.2 GB is the Rfam
CM database loading, and it is the same for either mode. It should be
set explicitly in the process `memory` directive.

**Against §3's acceptance criteria:**

- tRNA recovery ≥ 90 %: **not met** (78.4 % strict, 89.2 % lenient).
- < 10 % spurious: **met** (0 %).
- Both IR rRNA copies called: **met for 16S/23S/5S, not for 4.5S**.

§2 therefore **does not pass as written**. The shortfall is narrow and
predictable: every miss is an intron-containing tRNA or rrn4.5. That
is a fixed and well-known set of 8 + 2 loci in angiosperm plastomes,
and the same "intron structure" gap §2 already expected to leave open.
On everything else it is perfect.

### 7.2 Measurement — intron/*trans*-splicing in the existing miniprot baseline

Checked `INT-PLANT-01-pt/organelle_annotation.gff` against the lifted
`NC_018117.1` exon coordinates. The gap is smaller than §0 assumed:

- **Correct multi-exon structure:** `clpP` (3 exons), `ycf3` (3),
  `atpF`, `rpoC1`, `ndhB` (both IR copies), `rpl2`, `rps16`. Exon
  boundaries agree within a few codons.
- **Micro-exons missed** (protein alignment cannot place them):
  `petB` (6 bp exon 1), `petD` (9 bp), `rpl16` (9 bp). Each is called
  as its single large exon.
- **Spurious extra exon:** `ndhA` gets a third exon (120331–120436)
  that the reference does not have.
- **`rps12` *trans*-splicing:** only the IR-resident 3′ exons are
  called (protein residues 42–115, both IR copies). The 5′ exon in the
  LSC (71651–71762) is not called. This is as §0 predicted.

Neither §2 option touches this. Chloë is the only candidate that would.

### 7.3 Chloë (§4 items 1–2) — not yet asked

Re-checked 2026-10-13: neither `ian-small/Chloe.jl` nor
`ian-small/chloe_references` has a licence (GitHub API `license: null`,
`LICENSE` and `LICENSE.md` both 404). Last push dates are unchanged.
The email has **not** been sent; it needs a named DAFF sender. Draft:

> **Subject:** Licence for Chloë / chloe_references — use in a DAFF
> biosecurity pipeline
>
> Dear Prof. Small and Mr Castleden,
>
> We are building an open-source organelle-genome assembly and
> annotation workflow for the Australian Department of Agriculture,
> Fisheries and Forestry, used to identify intercepted plant and
> animal material from Nanopore reads. Chloë is the best fit we have
> found for plastome annotation, but neither the Chloe.jl nor the
> chloe_references repository carries a licence, so we cannot
> redistribute it in a container image.
>
> 1. Would you be willing to add an open-source licence (e.g. MIT or
>    Apache-2.0) to both repositories, or alternatively grant DAFF
>    written permission to redistribute them in a container image?
> 2. The README notes Chloë is optimised for angiosperms. How does it
>    behave on gymnosperm (e.g. conifer) and fern plastomes? Our
>    submissions are taxonomically unknown, and imported timber and
>    conifer nursery stock are routine.
>
> Many thanks,
> [name, role, DAFF]

### 7.4 Decision (2026-10-13): option A, `plant_pt` only

Decided by the user: **adopt barrnap alone for `plant_pt`**. `plant_mt`
stays CDS-only (shelved task 37). Implemented by task
65_plastid_noncds_barrnap.md, which also writes the spec decision
(§5's exit criteria for spec §3.1/§3.2 and spec §8 item 3). Chloë stays
a possible later upgrade. The §7.3 email is still unsent, which is the
only §5 criterion left open.

Options considered, given §7.1:

- **A. Adopt barrnap 1.10.6 alone** (`--kingdom bac --rrna --trna`).
  One bioconda container, no `mulled-build`, 6 s. Delivers 29 tRNAs +
  6 rRNAs, with zero spurious calls. The summary would record the
  known gaps: 8 intron-containing tRNAs, rrn4.5, rps12 5′ exon,
  micro-exons. This is a deviation from §3's threshold and from §2's
  two-tool bundle. It serves task 37 (`plant_mt`) with the same image.
- **B. tRNAscan-SE + barrnap as briefed.** Same strict recovery as A.
  The only extra is 4 low-score partial calls on intron tRNAs, which
  have the wrong structure and would need filtering out anyway. That
  means a `mulled-build` and about 14× the wall time for no net gain.
  Dominated by A.
- **C. `ok_cds_only` stands** pending the Chloë reply and a client
  answer to §4 item 3.

### 7.5 Other findings (outside this task's scope)

Briefed 2026-10-13: `rpl2` and the low scores in task
63_cds_locus_clustering.md, the mislabelled panel records in task
64_protein_panel_label_check.md.

- **`rpl2` IR copy B is dropped between `MINIPROT_CDS` and the final
  GFF.** `cds.gff` has good `rpl2` hits at 129890–131376; the merged
  `organelle_annotation.gff` keeps only the IR-A copy. Every other IR
  gene keeps both copies. One raw hit (`MP000542`, `NC_062422.1_rpl2`)
  spans 85769–131376, a 45 kb "spliced" alignment, and is a plausible
  culprit.
- **`clpP` is double-called as `ycf3`.** The panel entry
  `NC_034951.1_ycf3` is 196 aa (the length of clpP; ycf3 is ~168) and
  aligns exon-for-exon to the `clpP` locus. That looks like a
  mislabelled panel sequence, which puts a phantom second `ycf3` in
  the plastid annotation.
- **Low `qcovhsp` confidence on several correctly-spliced genes**
  (`atpF` 23, `rpoC1` 28, `ndhB` 27, `clpP` 28) while other multi-exon
  genes score ~100 (`ycf3`, `rps16`). Cause not investigated.
- §6's link to "task 33" was stale and has been corrected to
  task 37 (37_annotate_plant_mt_noncds.md, in `tasks/shelved/archived/`).

## 8. Scope addition (2026-10-13): plastid map visualisation

Added by the user: is there a UI library or tool for viewing an
annotated plastid genome, or should the existing organelle map
(task 44) be extended to cover plastids?

### 8.1 What the existing map already does for `plant_pt`

`bin/organelle_map.py` already renders `plant_pt`. On `INT-PLANT-01-pt`
it draws the single circular target contig (`path1`) plus the `path2`
SSC-inverted isoform, as inline SVG with a native `<title>` tooltip per
feature. "Display only the circular target contig" is therefore already
the behaviour. The gaps are in what a plastid reader expects to see:

| Gap | Effort | Notes |
|---|---|---|
| tRNA/rRNA features | **none** | `parse_non_cds_features` accepts any `tRNA`/`rRNA` GFF row with `Name=`. barrnap output would appear once ANNOTATE merges it. Its names (`tRNA-Met`, `16S_rRNA`) should be normalised to `trnM-CAU` / `rrn16` at merge time. |
| LSC / IRb / SSC / IRa segment ring | small | Segment lengths are already in `bin_metadata.json` (`lsc_len`/`ir_len`/`ssc_len`, task 44). This is the defining visual of a plastome map and the cheapest win. |
| Exon structure | small | Each multi-exon gene is drawn as one wedge across its mRNA span, introns included. The CDS children are already in the GFF, so draw exon blocks plus a thin intron connector. |
| Label every gene, not only barcodes | **not needed** | Decided by the user (2026-10-13): non-barcode genes are identified by tooltip only, as on `animal_mt`. This is already the behaviour: every feature on the `INT-PLANT-01-pt` map carries a `<title>` tooltip (161 on the fixture), and printed labels stay limited to barcode loci. |
| Zoom | hard in static SVG | A 72 bp tRNA is 0.17° of a 155 kb circle, about 1 px wide at report scale: visible, but fiddly to hover now that tooltips are how non-barcode genes are identified. Static SVG cannot fix this without JS. |
| GC-content ring, functional-group colours | small each | Needs `target.fasta` as a new input, and a gene→group asset. Nice-to-have. |

### 8.2 Off-the-shelf options

| Tool | Form | Verdict |
|---|---|---|
| **OGDRAW** (CHLOROBOX, MPI-MP) | Web server. A Perl module of the core was released in 2007. GenBank/EMBL input only. | The de facto standard look, but web-only for current versions ([rule 11](../CONSTITUTION.md)), GenBank input (we emit GFF3, spec §7), static image. **Rejected.** |
| **Chloroplot** (Zheng et al. 2020) | R/Shiny web app, source downloadable. GenBank input. | Plastid-specific (IR/SSC, codon usage), but static image, GenBank input, R stack. **Rejected** on task 44's own criteria. |
| **CPGview** | Web service. | Draws cis/*trans*-splicing maps (`rps12`), which is the most plastid-specific of the three, but web-only and static. **Rejected.** |
| **CPJSdraw**, IRscope | Perl / Shiny | IR-junction plots only, not a genome map. Out of scope. |
| **pyCirclize, plotgenes, Circos** | Python / R / Perl | Already rejected in task 44 (static image, no DOM). |
| **CGView.js** (Stothard lab; engine behind Proksee) | JS library, npm `cgview` 1.8.2 (2026-05-04), **Apache-2.0**, JOSS paper. Deps: `d3` 7, `svgcanvas`. `cgview.min.js` 222 KB plus d3. | **The one real candidate.** Interactive circular *and* linear views, zoom and pan to base level, built-in label collision handling, feature tracks/legends, GC tracks computed from the sequence, SVG export. Input is a JSON map that Python could produce directly from our GFF3 + `bin_metadata.json` (segments as a feature track). Inlining it keeps the report a single offline file, at about 0.5 MB per report (Plotly is already inlined). |

All the plastid-specific tools are web services or static-image
generators with GenBank input. They fail the same three tests that
ruled out pyCirclize et al. in task 44: no CLI, no DOM-level tooltips,
wrong input format. No offline, embeddable, plastid-specific map
library exists.

### 8.3 Options

- **V1. Extend the in-house SVG** (recommended first step). Add the
  segment ring, exon blocks, and tRNA/rRNA name normalisation. All are
  small, use data already in hand, and benefit `animal_mt` too (exon
  blocks are no-ops there). Zoom is left aside. Hover tooltips
  already name every feature.
- **V2. Adopt CGView.js for the report's interactive map.** This is
  the only route to zoom, which keeps small features (tRNAs) easy
  to hover. It is a real dependency decision ([rule 19](../CONSTITUTION.md)), and is
  best made for all three targets at once rather than for plastids
  only. It reverses task 44's "no JavaScript" choice for this one
  panel. The standalone `<id>.map.svg` could stay in-house.
- **V3. Both:** V1 now, V2 later if operators find the ~37 tRNA
  arcs too small to hover in practice.

**Decision (2026-10-13): V1.** Broken out as task
66_plastid_map_structure.md (segment ring, exon blocks, and
exon-overlap clustering). barrnap name normalisation stays with
whichever task implements §7.4's annotator decision. V2
(CGView.js) is deferred.

## Outcome (2026-10-15)

The user dropped the Chloë licence email (§7.3). No enquiry will be sent; Chloë stays recorded as a possible upgrade only if a licence appears. Task 65 implemented §7.4 (barrnap, `plant_pt`).
