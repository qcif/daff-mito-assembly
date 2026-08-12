# Task 34 — One-off polish-depth experiment on `CLIENT-BC05`

**Phase:** Out-of-phase. A **measurement task**, not a build task — it
changes no workflow code and ships no pipeline behaviour.
**Answers:** [spec §7 open question 6](../spec/07-open-questions.md)
("Polish depth (Flye built-in ± medaka)"), partially.
**Decides:** whether `tasks/shelved/17_medaka.md` gets un-shelved.

---

## Overview

Nanopore reads get most bases right but systematically miscount runs of
the same base — a real `AAAAAA` is read as `AAAAA` or `AAAAAAA`. In a
protein-coding gene a single inserted or deleted base shifts the reading
frame, so every downstream codon is garbage until another indel shifts it
back. "Polishing" is the corrective step: revisit the draft assembly with
the raw reads and let the consensus of many reads vote each base back to
the truth.

We currently polish exactly once, with Flye's own built-in polisher, and
we deliberately do not run `medaka` — the dedicated ONT polisher — because
`medaka` costs real wall-clock time and modern Dorado SUP reads are
already high-quality enough that it may buy nothing. That "may" has never
been measured on our own data. [Spec §7 item 6](../spec/07-open-questions.md)
records the open question; `MEDAKA` sits deferred at
[spec §2 stage 8](../spec/02-stages.md) with a numbered-but-empty slot
and a shelved implementation task, waiting for exactly this evidence.

A client sample has now given us a good place to measure. `CLIENT-BC05`
assembled into a single circular 15,908 bp animal mitogenome that BLASTs
at **99.08% nucleotide identity to *Bactrocera dorsalis*** (`NC_008748.1`,
Oriental fruit fly). That near-conspecific RefSeq record is the thing that
makes this experiment worth running: it is a curated, independent answer
key. We can align each polished variant against it and simply count how
many indels remain. That measures a polisher against a known answer rather
than tuning a knob until a fixture goes green — the distinction
[spec §5.1](../spec/05-test-data.md) insists on, and the same justification
[spec §7 item 12](../spec/07-open-questions.md) already accepts for
annotation benchmarking.

**Set expectations honestly before starting: the most likely outcome of
this experiment is "no meaningful gain".** The draft is already in good
shape — all 12 recovered protein-coding genes translate without a single
internal stop codon under genetic code 5, and all 6 barcode loci pass ORF
validation. [Spec §7 item 6](../spec/07-open-questions.md) scopes its
benchmark to "samples that fail ORF validation", and this sample does not
fail. A null result is therefore a *likely and completely valid* outcome,
and it is worth recording: it converts "we never checked" into "we checked
on modern SUP data and Flye's free pass was sufficient", which is the
evidence needed to keep `MEDAKA` deferred with confidence rather than by
assumption.

There is a second reason to run it anyway. Three narrow annotation
blemishes survive on this assembly (§2), all of the kind an indel would
explain. If one polish pass clears them, that is a concrete, mechanistic
win worth knowing about even if the aggregate ORF pass rate does not
move.

### Why this data is worth the effort

The reads self-report `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` in their FASTQ
headers, from a run dated 2026-04-15. That matters beyond this experiment:
the [`tasks/todo.md`](todo.md) "Benchmark data" item records that we have
**no** modern benchmark set — the Tier 2 integration fixtures are old
public SRA predating the R10.4.1 SUP chemistry the pipeline assumes. This
client sample *is* R10.4.1 SUP at realistic depth. It does not close that
gap on its own (one taxon, one target, `animal_mt` only), but it is the
first data we have that matches the chemistry the pipeline was designed
around, and any figure measured here is more representative than a fixture
figure.

---

## 1. Constraints

**No workflow changes. None.** This task does not touch `main.nf`,
`modules/`, `conf/`, `nextflow.config`, `bin/`, `assets/`, or
`tests/`. `params.polish` stays `false`;
[`modules/local/medaka.nf`](../modules/local/medaka.nf) keeps its stub.
The experiment runs as standalone containerised commands against files
that already exist on disk. The only committed artefacts are this file
(with §7 filled in) and any spec/todo updates §8 calls for.

Rules from [CONSTITUTION.md](../CONSTITUTION.md) that still bind, because
they bind everywhere:

- **Rule 11 — every tool runs in a pinned container.** A one-off is not a
  licence for host tools or a conda env. Flye and medaka both come from
  SHA-pinned biocontainers, and the digests get recorded in §7. A number
  produced by an unpinned tool is not auditable and does not count as
  evidence.
- **Rule 12 — prefer biocontainers.** Both tools have published ones; no
  bespoke image, no `mulled-build`.
- **Rule 18 — favour auditability.** Every arm records its exact command,
  container digest, input file, and wall time. An unreproducible
  measurement cannot justify un-deferring a stage.
- **Rule 19 — maintenance.** The deliverable is a written finding, not
  new machinery. Do not build a benchmark harness for a single sample.

Per the project's Docker convention, pass `--user $(id -u):$(id -g)` on
every interactive `docker run` so outputs are not written as root.

**Not in scope:**

- Un-deferring `MEDAKA`, editing `tasks/shelved/17_medaka.md`, or setting
  `params.polish = true`. This task produces the evidence that decides
  those; it does not act on it. Acting on it is a follow-up task, and if
  the finding is positive that follow-up is the already-written shelved
  task 17, unshelved and renumbered.
- Iterative polishing (2–3 rounds of either tool). Spec §7 item 6 lists
  iterative arms (b) and (d); this task runs **one pass only** per arm, by
  explicit instruction. If single-pass shows a real gain, the iterative
  arms become a follow-up; if it shows nothing, they are moot.
- The other two assembly targets. `plant_pt` and `plant_mt` have no
  comparable near-conspecific answer key in hand and no client sample run
  yet.
- Re-running the pipeline, re-assembling, or altering `CLIENT-BC05`'s
  published outputs under `tests/manual/output/`. Arm A is read-only.
- Generalising from n=1. See §6.

---

## 2. Baseline — what we are actually trying to improve

From the completed `CLIENT-BC05` run
(`tests/manual/output/CLIENT-BC05/`). Record these numbers in §7 as the
Arm A column; re-derive them rather than trusting this table, since the
run predates the task.

| Signal | Baseline value |
|---|---|
| Target contig | `contig_1`, 15,908 bp, circular (`flye_circ`), 65× |
| Protein-coding completeness | 0.923 — 12/13, **`ATP8` missing** |
| tRNA / rRNA | 21 tRNA (of 22 — `trnS2` absent), 2 rRNA |
| Internal stop codons, all 12 CDS, table 5 | **0** |
| Barcode ORF validation | **6/6 pass** (COX1, COX2, COX3, CYTB, ND1, ATP6) |
| CDS carrying a miniprot `Frameshift` attribute | **1 of 12 — `ATP6`, `Frameshift=2`** |
| miniprot/MITOS2 coordinate conflicts | **1 — `nad4`** |
| Best-hit protein identity range | 0.54–0.80 |

**Three blemishes are the real targets of this experiment**, and they are
the only things a polish pass could plausibly fix:

1. **`ATP6` carries `Frameshift=2`.** The single unambiguous piece of
   indel evidence in the annotation. If one polish pass drops this to 0,
   that is a direct, mechanistic demonstration that polishing repairs
   real frame damage on our data.
2. **`ATP8` is missing entirely.** `ATP8` is short and fast-evolving and
   is commonly missed by miniprot for ordinary reasons of divergence, so
   this is *probably* panel distance rather than damage — but a frameshift
   early in a short gene is also a sufficient explanation, and polishing
   discriminates between the two.
3. **`nad4` coordinates disagree between miniprot and MITOS2.** An indel
   near a gene boundary shifts where each method calls the start/stop, so
   an annotation conflict is a plausible indel signature.

**Do not treat low protein identity (0.54–0.80) as damage.** It is
expected panel distance, already diagnosed in
[task 30 §10.1](completed/30_unified_locus_pass.md) and carried in
[`tasks/todo.md`](todo.md): the `animal_mt` protein panel spans all of
Metazoa at 10 representatives per gene and contains no tephritid fruit
fly. A 99.08% nucleotide identity to `NC_008748.1` alongside 54–80%
protein identity to a distant panel is entirely self-consistent and is a
panel-breadth artifact, not sequencing error. Polishing will not move it
much, and if §7 reports a large identity jump, suspect the measurement
before believing it.

> **Correction of record.** An earlier informal reading of this run
> interpreted miniprot's `StopCodon=N` GFF attribute (values 1–14 across
> the CDS) as premature in-frame stops, and concluded the assembly was
> substantially indel-damaged. **That was wrong.** Direct translation of
> all recovered CDS under table 5 yields **zero** internal stops, and
> `validate_barcodes.py` passes 6/6. Whatever miniprot's `StopCodon`
> attribute counts, it is not premature stops in our sequence. Confirm
> miniprot's actual documented semantics for that attribute while running
> this task and note it in §7 — it appears in every CDS line we emit and
> we should not misread it a second time.

---

## 3. Experimental design

Three arms, one pass each, same contig and same reads throughout.

| Arm | What | Cost |
|---|---|---|
| **A — baseline** | `contig_1` exactly as it stands. Flye's built-in polish at its `--nano-hq` default of `--iterations 1`, which is what [`modules/local/metaflye.nf`](../modules/local/metaflye.nf) produces (it passes no `--iterations`, so Flye's default applies). | Already computed — read from disk, do not re-run. |
| **B — one extra Flye pass** | Flye's standalone polish mode (`--polish-target`) applied to Arm A's contig with `--iterations 1`. Isolates "is Flye's own polisher simply under-run?" — the free option, and per spec §7 item 6 the baseline any paid option must beat. | Minutes. |
| **C — one medaka pass** | `medaka_consensus` applied to Arm A's contig. The paid option. This is exactly what deferred stage 8 would do, so a null result here is a direct answer about stage 8. | Tens of minutes, CPU. |

Arms B and C both start from Arm A — they are alternative *additions* to
the shipped pipeline, not alternatives to each other stacked in sequence.
Do not feed B's output into C.

### Inputs — use the reads METAFLYE actually consumed

Both polishers must see the same evidence the assembler saw, or the
comparison is meaningless. That is the **coverage-gated** FASTQ, which is
`COVERAGE_GATE`'s output and `METAFLYE`'s input per
[main.nf](../main.nf) — not the raw sample and not the recruited pool.
For this sample the gate reported `status: "ok"` at 143.35× against a
`max_cov` of 300, so it passed reads through without subsampling; the
gated FASTQ is therefore equivalent in content to
`tests/manual/output/CLIENT-BC05/recruit/CLIENT-BC05.recruited.fastq.gz`
(576 reads). Confirm that equivalence rather than assuming it — if the
`METAFLYE` work directory still exists, take its staged input file
directly and prefer it.

The contig to polish is the **binned target**, not the whole assembly:
`tests/manual/output/CLIENT-BC05/bin_target/target.fasta` (`contig_1`
alone). Polishing all 8 contigs would waste time on off-target material
that `BIN_TARGET` already classified away.

### Medaka model

No guesswork needed, and this resolves the open question shelved task 17
§1 flagged. The reads carry
`basecall_model_version_id=dna_r10.4.1_e8.2_400bps_sup@v5.0.0` in their
FASTQ headers, which maps to medaka model **`r1041_e82_400bps_sup_v5.0.0`**
— the exact string shelved task 17 proposed as the default. Verify it is
present in `medaka tools list_models` inside the pinned image before
running; if 2.2.2's catalog has dropped it, take the nearest R10.4.1 SUP
v5 model and record the substitution and reason in §7.

### Containers

Pin both by digest and record both in §7:

- Flye — the digest already pinned for `METAFLYE` in
  [`conf/containers.config`](../conf/containers.config)
  (`flye:2.9.6--py311h93bbee8_1@sha256:eb57665b…`). Reusing the shipped
  pin means Arm B measures *our* Flye, not a different build.
- medaka — `quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0`, the tag
  shelved task 17 §2 nominates. Resolve its digest before first use and
  record it; a follow-up that un-shelves task 17 can then reuse the same
  digest rather than re-resolving.

Illustrative shape only — resolve real paths and digests at run time:

```bash
# Arm B — Flye standalone polish, one iteration
docker run --rm --user $(id -u):$(id -g) -v "$PWD":/w -w /w \
  <flye-image@sha256:...> \
  flye --polish-target <target.fasta> \
       --nano-hq <gated.fastq.gz> \
       --iterations 1 \
       --threads <n> --out-dir armB_out

# Arm C — medaka, one pass
docker run --rm --user $(id -u):$(id -g) -v "$PWD":/w -w /w \
  <medaka-image@sha256:...> \
  medaka_consensus -i <gated.fastq.gz> -d <target.fasta> \
       -m r1041_e82_400bps_sup_v5.0.0 -o armC_out -t <n>
```

Work under a scratch directory outside the repo, or a gitignored path.
Do not write into `tests/manual/output/`, which holds Arm A.

---

## 4. Metrics

### 4.1 Primary — indels against the answer key

The measurement this experiment exists for. Align each arm's contig to
**`NC_008748.1`** (*Bactrocera dorsalis*, the 99.08% BLAST hit) and count
indel events and total indel bases.

Fetch `NC_008748.1` from NCBI and record its accession.version plus
retrieval date in §7 — it is external reference data and per
[rule 10](../CONSTITUTION.md) an untracked reference makes the result
unauditable. It is a one-off download for a one-off experiment and does
**not** enter `refs/` or `manifest.json`.

Align with a pinned aligner already in the project's container set
(`minimap2` at `asm5`/`asm10` suits two near-identical mitogenomes) and
count indel operations from the CIGAR. Both contig and reference are
circular with arbitrary start points, so **rotate or double the reference
before aligning** or the alignment will artificially break at the origin
and manufacture apparent indels. State in §7 how origin was handled.

Report per arm: indel events, total indel bases, and how many fall inside
homopolymer runs of ≥4 — the last being the specific error mode polishing
is supposed to fix, and therefore the most diagnostic number in the whole
experiment.

**Interpretation caveat, to be stated in §7 rather than glossed:**
99.08% identity means ~0.92% genuine divergence, so not every difference
from `NC_008748.1` is our error — *B. dorsalis* is a close relative, not
the same individual. Absolute indel counts therefore overstate assembly
error. The comparison across arms on the *same* reference remains valid,
because real biological divergence is constant across arms and cancels;
only the *differences* between A, B and C are attributable to polishing.

### 4.2 Secondary — does the annotation improve?

Re-run the annotation and barcode logic on each arm's contig and diff
against Arm A. Do this by invoking the same containerised tools the
pipeline uses, against the arm's FASTA — **not** by re-running the
workflow.

- `ATP6` `Frameshift` attribute: 2 → 0? *(the sharpest single signal)*
- `ATP8`: recovered, taking completeness 12/13 → 13/13?
- `nad4` miniprot/MITOS2 coordinate conflict: resolved?
- Internal stop codons across all CDS under table 5: still 0? **A polish
  pass that introduces a stop is a red flag, not a neutral result** —
  report it prominently if it happens.
- Barcode ORF validation: still 6/6?
- Per-gene best-hit protein identity: direction and magnitude of change.
  Expect small movement (see §2's warning); a large jump means measurement
  error.

### 4.3 Cost and sanity

- Wall-clock and peak RSS per arm, on stated CPU/RAM. Spec §7 item 6's
  success metric is explicitly "ORF pass rate **vs. wall-clock cost**" —
  a gain that costs 40 minutes on a 16 kb mitogenome is a different
  proposition from one that costs 4.
- Total polished length within **±5%** of Arm A (shelved task 17 §4's
  sanity check). Polishing edits homopolymers; it must not add, drop or
  fragment contigs. A violation means the run is wrong, not that the
  polisher is bad.
- Contig count stays 1 and circularity is not silently lost.

---

## 5. Decision rule

Fix this before looking at the numbers, so the conclusion is not fitted to
the result.

- **Un-defer `MEDAKA`** (i.e. unshelve task 17) only if Arm C beats both A
  and B on §4.1 indels *and* delivers at least one §4.2 improvement, at a
  wall-clock cost proportionate to a 16 kb assembly. Per spec §7 item 6,
  medaka must beat **Flye's free extra pass**, not merely the baseline —
  if B and C land in the same place, B wins on cost and dependency count
  ([rule 19](../CONSTITUTION.md)).
- **Consider raising Flye `--iterations`** if B clearly beats A and
  roughly matches C. That is a one-line change to
  [`modules/local/metaflye.nf`](../modules/local/metaflye.nf) with no new
  dependency — the cheapest possible win. It would still need its own task
  and its own spec §7 item 6 note; do not make the edit here.
- **Keep the status quo** if neither arm moves §4.1 or §4.2 materially.
  This is the expected outcome (§Overview) and is a **successful result**,
  not a wasted run. Record it and move on.
- **Escalate** if any arm makes things worse — a new internal stop, lost
  circularity, a length violation. That would mean polishing is actively
  harmful on this data and belongs in spec §7 as a finding in its own
  right.

---

## 6. Limits on what this can conclude

State these plainly in §7 so the finding is not over-read later.

- **n = 1.** One sample, one taxon, one target, one contig. It cannot
  establish a general polish policy — it can only justify keeping stage 8
  deferred, or motivate a larger benchmark.
- **This sample is near the quality ceiling.** 6/6 ORF pass, 0 internal
  stops, 65× coverage, R10.4.1 SUP reads. Spec §7 item 6 asks for
  "samples that fail ORF validation" and this is not one. A null result
  here says *"polish adds nothing on good data"* — it does **not** say
  polish adds nothing on marginal data, which is precisely the case
  `--polish` exists to serve. Say so explicitly.
- **The answer key is a relative, not a conspecific** (§4.1 caveat).
- **`animal_mt` only.** Plant arms are untested and plant mitogenomes are
  larger, repeat-rich and structurally harder; nothing here transfers to
  them.

---

## 7. Outcomes

_To be filled in on completion._ Record, at minimum:

- Container digests actually used (Flye, medaka, aligner), and the medaka
  model string resolved.
- `NC_008748.1` accession.version + retrieval date; how circular origin
  was handled in §4.1.
- The §4.1 / §4.2 / §4.3 table, three arms side by side.
- Miniprot's documented semantics for the `StopCodon=` GFF attribute
  (see §2's correction of record).
- The §5 decision reached, and the §6 limits restated.

---

## 8. Follow-through

Once §7 is filled in, whichever way it lands:

- Update [spec §7 item 6](../spec/07-open-questions.md) — replace the
  unqualified "MEDAKA stage deferred" framing with the measured result and
  a pointer to this task ("task 34_polish_depth_oneoff.md"). If the
  decision is to keep deferring, say *measured on modern SUP data*, which
  is stronger than the current assumption-based deferral.
- Add a line to [`tasks/todo.md`](todo.md) under **Benchmark data**
  noting that `CLIENT-BC05` is the first R10.4.1 SUP sample we have, that
  it is `animal_mt`-only, and what it can and cannot support. That item
  currently states we have no modern data at all; this partially amends
  it.
- If the decision is to un-shelve, move `tasks/shelved/17_medaka.md` into
  `tasks/` and renumber it to the next free number, folding this task's
  §7 numbers into its brief as the justification. Its §1 medaka-model
  question is already answered above.

## 9. Test surfaces

**None, deliberately.** [Rule 15](../CONSTITUTION.md) defines three test
surfaces for workflow code; this task ships no workflow code, no
`bin/*.py`, and no module change, so none of the three apply. In
particular:

- No `pytest` — [rule 14](../CONSTITUTION.md) governs custom logic
  promoted into `bin/`, and nothing here is promoted. Any scratch
  scripting for counting indels or diffing GFFs is throwaway analysis
  under a scratch path, not project code. **If that scratch work grows
  past ~10 lines and looks reusable, that is a signal to stop and write a
  proper task**, not to quietly add an untested script to `bin/`.
- No integration fixtures. No stage behaviour changes, so
  `tests/integration/expected/*` and `assertions.sh` stay untouched, and
  no new fixture is staged — `CLIENT-BC05` is client data and is **not**
  a candidate for the public fixture set without a separate provenance and
  data-sharing review.
- `-stub-run` and `-profile integration` must remain green *by virtue of
  nothing having changed*. If either breaks, something was edited that
  §1 forbids.
