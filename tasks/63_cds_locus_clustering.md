# Task 63 — CDS locus clustering: exon overlap, strand, and an intron cap

**Phase:** P4b (plant annotation arms). Rolled out of task
32_research_plastid_noncds.md §7.2 and §7.5.

**Status:** in progress on branch origin/worktree-task-63. Last update:

>Task 63 is partly done: the code, unit tests and spec updates are committed, but the integration side needs a machine with Nextflow. I committed to the task-63 worktree (a40091a, branch worktree-task-63). I did not merge it, and the task file is still in tasks/ rather than tasks/completed/.
>Done
>
>Shared clustering: cluster_loci in bin/annotation_gff.py now groups hits by gene, contig and strand, and links them only when an exon overlaps an exon. Both the barcode extractor and the annotator call it, so they pick the same winner. Ties go to the leftmost start, then the lowest mRNA ID.
>Intron cap: params.miniprot_max_intron is 5,000 bp for plant_pt and 200,000 for the other two targets. It is passed to miniprot as -G and recorded as max_intron in genetic_code.json. The animal trial step (select_genetic_code.py) now needs a --max-intron argument to write it.
>Tests: I added unit tests for the clustering, including the rpl2 case and a check that the barcode and annotation callers agree. 293 tests in the touched and neighbouring test modules pass.
>Integration assertion: plant_pt/annotation_bounds.json now lists ir_pair_genes, and assertions.sh checks each appears at least twice. I only checked that the script and JSON are syntactically valid.
>Spec and todo: I updated the stage 12, barcode extractor and annotator rows, added a new §9 item 13 for the intron cap, and added the barcode length check to tasks/todo.md.
>Not done, because this machine has no Nextflow, Docker or fixture outputs
>
>The stub run (§5.2).
>The before/after integration comparison (§5.3), which is where the expected rpl2 1 → 2 change would show up.
>The frameshift diagnosis (§2.3). So the spec §9 item 6 text and the possible Annotation todo item are still to write.
>The 100% branch coverage and flake8 checks, since neither tool is installed here.
>I reported this wrongly earlier. At one point I said the new unit tests were written, but that edit had been interrupted. They are written now and run. The "Not done" items above are also listed in the task's Outcomes section.
>Decision for you: run the outstanding checks in CI or on the workstation before merging, or merge now and track them separately? tasks/active.md on main still lists task 63. It's uncommitted, and it should stay until the merge.


**No preconditions.** Independent of task
64_protein_panel_label_check.md and task 65_plastid_noncds_barrnap.md.
Should land before task 66_plastid_map_structure.md, whose §3.3 makes
the same clustering change in the map and should reuse this task's
predicate rather than write a third copy.

## Overview

Stage 12 (`MINIPROT_CDS`) aligns every protein in the reference panel
against the assembled organelle. A real gene usually attracts several
alignments, one per panel representative (about ten per gene for
`plant_pt`). Two custom components then collapse those into one call
per gene locus:

- C5 (`bin/validate_barcodes.py`, `cluster_by_locus`) for the barcode
  FASTA;
- C8 (`bin/annotate_summary.py`, `cluster_cds_by_gene`) for the
  published annotation GFF.

Both group hits of the same gene name whose **outer span** (first exon
start to last exon end) overlaps, ignoring strand, and keep the
highest-identity hit in each group. A gene that appears at two places
in the genome stays as two features, which is how a plastid keeps both
copies of its inverted-repeat (IR) genes.

On the plastid fixture this goes wrong. miniprot produced one
nonsensical `rpl2` alignment whose "intron" is 44 kb long. Its outer
span reaches from one IR copy of `rpl2` to the other. Span clustering
chains every `rpl2` hit on either strand into one group, and the
published annotation keeps one `rpl2` where the genome has two.

This task fixes the problem in two places:

1. Hits group only when they are on the **same strand** and their
   **exons** overlap. A hit sitting inside another hit's intron no
   longer joins its group.
2. miniprot gets a **per-target maximum intron size**, so an
   alignment like the 44 kb `rpl2` is never produced for `plant_pt`.

There is one shared clustering routine, not two. C5 and C8 have to
agree on what a locus is: the barcode cut from a locus must be the
same feature the annotation shows (task 30's barcode-coherence
invariant). Today they agree only because the code was copied.

The task also settles task 32's third finding: **low confidence
scores** on correctly spliced genes. It turns out not to be about
splicing (§1.3). Here that finding is diagnosed and recorded, not
fixed.

## 1. Evidence (`INT-PLANT-01-pt`, 2026-10-13)

From `tests/integration/output/INT-PLANT-01-pt/` on the current
bundle (`v2026.09_1`).

### 1.1 The lost `rpl2` copy

Raw `cds/INT-PLANT-01-pt.cds.gff` has two clean groups of `rpl2`
hits, one per IR copy, plus one bridging hit:

| Hits | Strand | Span | Exons |
|---|---|---|---|
| 9 (MP000531 … MP000548) | + | 109,882–111,380 | 109,894–110,283 / 110,937–111,380 |
| 9 (MP000530 … MP000547) | − | 129,890–131,388 | 129,890–130,333 / 130,987–131,376 |
| **MP000542** | − | **85,769–131,376** | 130,987–131,376 / 129,901–130,333 / **85,769–85,881** |

MP000542's last "exon" is 113 bp at 85.8 kb, joined by a 44 kb
intron. `cluster_cds_by_gene` sorts by start, opens a group at
85,769, and extends the group's end to 131,376. Every later `rpl2`
hit then falls inside the group, whatever its strand. The winner is
MP000546 (identity 0.958, + strand), and the − strand IR copy
disappears from `organelle_annotation.gff`. All other IR genes
(`ndhB`, `rps7`, `rpl23`, `ycf2`, `rps12`) keep both copies.

Exon-plus-strand clustering alone would fix this instance: MP000542
joins the − strand group, then loses on identity (0.9455 < 0.958). It
does not fix the general case. If a bridging hit ever has the best
identity, it wins and publishes a 45 kb gene. That is why §2.2 also
removes such hits at the source.

### 1.2 Plastid intron sizes

The real `plant_pt` multi-exon calls on the fixture have introns of
about 0.6–0.9 kb (`rpl2` 653 bp, `ycf3` ~0.7 kb, `clpP` ~0.6–0.8 kb,
`atpF` ~0.7 kb, `ndhB` ~0.6 kb). The longest cis-spliced intron in a
land-plant plastid protein-coding gene is about 1.1 kb (`ndhA`).
`rps12` is *trans*-spliced and miniprot never calls its 5′ exon (task
32 §7.2), so a cap cannot lose it. miniprot's default cap is 200 kb
(`-G`, confirmed in the pinned `miniprot:0.18--h577a1d6_0` image's
usage text). That is larger than the whole plastome.

### 1.3 Low confidence scores track frameshifts, not splicing

Task 32 §7.5 reported low `qcovhsp` on several correctly spliced genes
(`atpF` 23, `rpoC1` 28, `ndhB` 27, `clpP` 28). Across all 90 mRNA
features in the fixture's published GFF:

| miniprot `Frameshift=` attribute | Exons | Features | `qcovhsp` < 50 |
|---|---|---|---|
| absent | single | 31 | **0** |
| absent | multi | 1 | **0** |
| present | single | 48 | 30 |
| present | multi | 10 | 8 |

Every low score is on a feature where miniprot had to put a frameshift
in its alignment, and 48 of the 58 frameshifted features are
single-exon. The cause is the frameshifts, not splicing.
`translate_annotation_cds.py` translates the spliced sequence as one
reading frame by design (task 40 §1), so a frameshift scrambles
everything downstream of it. The low score is the intended result.

The frameshifts also reach the barcodes. All six `plant_pt` barcodes
pass ORF validation, and none has a length divisible by 3 (`rbcL`
1,451 nt against a ~1,428 nt reference CDS; `matK` 1,522; `rpoB`
3,205). ORF validation drops frameshift ops before translating
(`codon_blocks`), so it cannot see them.

The likely cause is uncorrected ONT indels, mostly in homopolymers.
`INT-PLANT-01-pt` is old public SRA data that predates the R10.4.1 SUP
chemistry the pipeline assumes (spec §5.1), and polishing is Flye's
single built-in pass (MEDAKA deferred, spec §9 item 6). One case
needs a closer look: `atpF`'s only frameshift (`41M2G9M…`) sits 9
codons before its splice junction, which can be a miniprot artifact
(junction misplaced by 2 nt, then compensated) rather than an
assembly error.

## 2. Scope

### 2.1 One shared locus-clustering routine (C5 + C8)

- Move the clustering into one stdlib-only function in a shared
  `bin/` module (e.g. `bin/intervals.py`, which already plays this
  role for merged intervals, or `bin/annotation_gff.py`). It must stay
  stdlib-only: C8 runs inside the MITOS2 biocontainer, which has no
  biopython (task 31 §2). Both `cluster_by_locus` (C5) and
  `cluster_cds_by_gene` (C8) become thin callers.
- **Grouping key:** gene symbol (each caller keeps its current
  case handling: C5 upper-cases, C8 does not), `seqid` and
  **strand**.
- **Overlap test:** two hits are in the same locus when **any exon of
  one overlaps any exon of the other**. Clusters form transitively.
  An exon-sorted sweep is fine; the inputs are tens of hits per gene.
- **Winner:** unchanged. Highest `identity`, ties broken
  deterministically (state the rule and test it).
- **Never de-duplicate by gene name** (task 31 §4 point 6). The new
  rule only ever *splits* groups the old rule merged. It never merges
  groups the old rule kept apart.
- Task 66_plastid_map_structure.md §3.3 makes the same change to
  `organelle_map.cluster_loci()`. Write the exon-overlap predicate so
  that task can import it. Do not edit `organelle_map.py` here.

### 2.2 Per-target miniprot intron cap

- New param, a per-target map in `nextflow.config` next to
  `genetic_code_tables` (principle 9: versioned config, not
  hard-coded), passed to miniprot as `-G`:

  | Target | Value | Reason |
  |---|---|---|
  | `plant_pt` | 5,000 bp | ~4.5× the longest known cis intron (§1.2); room for non-angiosperm lineages without allowing a cross-IR bridge |
  | `plant_mt` | 200,000 (miniprot default, written out) | Plant mitochondrial cis-introns are larger and less well characterised; no fixture evidence to tune on |
  | `animal_mt` | 200,000 (miniprot default, written out) | No introns. Lowering it is harmless but untested, so leave it |

- Writing the defaults out explicitly records the value in config and
  provenance, so an auditor doesn't need to know miniprot's default.
- The process reads the value as a `val`, the same way
  `genetic_code_tables` is used today. No new staged file.
- Record the cap applied in `genetic_code.json` (MINIPROT_CDS's own
  provenance file, already published) or in
  `annotation_summary.json`, whichever already carries miniprot's
  settings. Pick one and say which in Outcomes (rule 16).

### 2.3 Frameshift diagnosis (measure and record, no code change)

1. For the four genes task 32 named (`atpF`, `rpoC1`, `ndhB`, `clpP`)
   and the six barcodes, locate each miniprot frameshift on the
   fixture (the `cg:Z:` CIGAR on the `##PAF` line in `cds.gff`). Align
   the surrounding ~30 bp to the same gene in `NC_018117.1`, the
   fixture's 99.60 %-identity RefSeq match (task 21 §9). Classify
   each frameshift:
   - **assembly indel:** the assembly lacks or adds a base the
     reference has, typically in a homopolymer;
   - **alignment artifact:** the assembly matches the reference and
     the frameshift is miniprot's placement, typically at a splice
     junction.
2. **If every frameshift is an assembly indel:** no change. The score
   is doing its job. Record the finding in Outcomes and in spec §9
   item 6 (polish depth) as fixture evidence that Flye's single polish
   pass leaves frameshifts in old-chemistry data. Do **not** re-tune
   anything on this fixture (spec §5.1).
3. **If any are alignment artifacts at splice junctions:** record
   them and add a `tasks/todo.md` item under Annotation. The fix
   (translating each CDS segment in the frame its GFF phase column
   gives) changes the deliberate design of task 40 §1. It needs its
   own brief, not a rider on this one.
4. Either way, add one `tasks/todo.md` line under Barcode extraction:
   ORF validation passes barcodes whose length is not a whole number
   of codons, because it drops frameshift ops before translating. Note
   whether a `length_nt % 3` warning in `validation.tsv` would be
   worth adding. Do not add it here.

### 2.4 Out of scope

- The phantom `ycf3` on the `clpP` locus and the phantom `cemA` on the
  `petA` locus. These are mislabelled panel records, handled in task
  64_protein_panel_label_check.md. This task's clustering is per gene
  name and cannot see them.
- Any runtime guard against two *different* gene names called on the
  same exons. Fix mislabels at the source (task 64), not by adding a
  second arbiter at runtime.
- `rps12` 5′ exon and the micro-exons of `petB`, `petD` and `rpl16`
  (task 32 §7.2). Protein-to-genome alignment cannot place them.
- Trimming miniprot spans at frameshifts (`tasks/todo.md`,
  Annotation, task 39 item (b)).
- Changing the confidence-score translation (see §2.3 point 3).

## 3. Constraints

- **Barcode coherence (task 30, spec §8 item 3).** C5 and C8 must
  pick the same winner for every locus. Sharing one function
  guarantees it. A test pins it (§5.1).
- **CDS coordinates are never modified** (task 31 §4 point 1). This
  task changes which hits are kept, never their coordinates.
- **Panel-competition check (spec §4.3).** The intron cap changes
  miniprot's output, and the clustering change can change a winner.
  Neither may lose or degrade a barcode locus on any fixture (§5.3).
- **Principle 8 / C8 always exits 0.** No new failure path.

## 4. Pseudocode

```
def exons_overlap(a, b):
    return any(s1 <= e2 and s2 <= e1
               for s1, e1 in a.exons for s2, e2 in b.exons)

def cluster_loci(records, key=gene_symbol):
    groups = group_by(records, (key(r), r.seqid, r.strand))
    winners = {}
    for (gene, seqid, strand), hits in groups:
        clusters = connected_components(hits, edge=exons_overlap)
        winners[gene] += [max(c, key=(identity, tiebreak)) for c in clusters]
    return winners

# MINIPROT_CDS script block
miniprot --gff -T $table -G $max_intron ...
```

## 5. Testing

### 5.1 Unit tests (CONSTITUTION rules 14, 15)

New module or new tests next to the shared function. Keep C5
(`scripts/tests/test_validate_barcodes.py`) and C8
(`scripts/tests/test_annotate_summary.py`) at 100 % branch coverage.
Inputs are small GFF strings written to `tmp_path`, with nothing
mocked below the module boundary (task 27).

- **The `rpl2` shape:** two compact groups on opposite strands, plus
  a bridging hit whose last exon is far upstream. Expect two winners.
  Also run it with the bridging hit at the highest identity. It wins
  its own strand's group only, and the other strand's copy survives.
- **Strand:** identical exons on opposite strands give two loci.
- **Nested in an intron:** a same-gene, same-strand hit inside another
  hit's intron, with no exon overlap, gives two loci.
- **Exon overlap still merges:** redundant panel hits with slightly
  different exon boundaries give one locus, highest identity wins.
- **Transitive merge:** A overlaps B and B overlaps C, but A does not
  overlap C. Expect one locus.
- **Tie-break is deterministic:** equal identity gives the same winner
  regardless of input order.
- **Never de-duplicate by gene name:** the existing IR-copy test in
  `test_annotate_summary.py` still passes unchanged.
- **C5/C8 agreement:** the same `cds.gff` through both callers gives
  identical winners (IDs and coordinates) for every barcode gene.

### 5.2 Stub run

The new param is a `val`, and no channel shapes change. Confirm that
`-stub-run` passes, and that CI's bare-params lint passes with the
param in the script block.

### 5.3 Integration and panel competition

Run `-profile integration` before and after the change. Record, per
fixture:

- **Barcodes (contractual):** every locus in `validation.tsv` still
  `pass`, with the same coordinates. A coordinate shift is
  acceptable only with improved identity and query coverage, and must
  be explained in Outcomes (spec §4.3 panel-competition rule).
- **Annotation:** CDS feature count per gene before and after. Every
  gene that *gains* a copy must be explained, because splitting
  clusters can expose a same-gene fragment that span clustering used
  to swallow. Expected on `INT-PLANT-01-pt`: `rpl2` 1 → 2, nothing
  else changes.
- `INT-ANIMAL-01`: no change at all is expected (no introns, cap
  unchanged).

Add `tests/integration/expected/plant_pt/annotation_bounds.json`, in
the shape of the `animal_mt` file. Its first assertion is
`ir_pair_genes`: genes that must appear **at least twice** in the
published GFF. Use `rpl2`, `ndhB`, `rps7`, `rpl23`, `ycf2`. IR
duplication is biology, not a tuned value, so it won't drift with the
reference data (rule 19). Wire it into the existing `INT-PLANT-01-pt`
annotation block in `tests/integration/assertions.sh`. Task 64 and
task 65 add to this file. Create it here so they extend it.

## 6. Spec updates

- spec/02-stages.md, stage 12 row: the per-target `-G` cap and its
  reason. C5 and C8 rows: locus clustering is on same-strand exon
  overlap, in one shared function.
- spec/07-open-questions.md §9 item 6 (polish depth): the §2.3
  frameshift finding as fixture evidence, worded as in §2.3 point 2.
- spec/07-open-questions.md §9: add the intron cap to the list of
  provisional per-target values to sweep against the benchmark set.

## 7. Checklist

- [ ] §2.1 shared clustering in one module; C5 and C8 call it
- [ ] §2.2 `-G` per-target param wired and recorded in provenance
- [ ] §2.3 frameshifts classified; todo items added
- [ ] §5.1 unit tests; 100 % branch on C5, C8 and the shared module;
      flake8 clean
- [ ] §5.2 stub run passes
- [ ] §5.3 before/after recorded; no barcode lost or degraded;
      `plant_pt/annotation_bounds.json` created and asserted
- [ ] §6 spec updated
- [ ] Outcomes section written
