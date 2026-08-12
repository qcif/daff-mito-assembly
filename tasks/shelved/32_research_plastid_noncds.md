# Task 32 — Research: plastid (`plant_pt`) non-CDS annotation

**Phase:** P4 (from [spec §6](../spec/06-phases.md)), **research only —
no implementation.** Produces a decision and its evidence, not code.

**Status:** open, and **no longer blocking**. See §0.

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
- **`plant_mt` non-CDS is [task 33](33_annotate_plant_mt_noncds.md)**,
  blocked on its own preconditions. If §2 is accepted here, task 33
  becomes nearly free — the same two tools serve both plant targets.
