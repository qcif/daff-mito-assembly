# Task 44 — Stage 14 `ORGANELLE_MAP`: annotated organelle diagram

**Phase:** P4a (from spec §6).
**Goal:** Replace the P0 stub in
`modules/local/organelle_map.nf`
(currently `touch ${meta.sample_id}.map.svg`) with a real renderer
producing inline-ready SVG from `ANNOTATION_SCORING`'s GFF3 +
`annotation_summary.json`.

**Carries an undecided design question.**
spec §8 item 4
marks the tool choice deferred and says the decision is "needed before
P4a" — so **this task opens by making that decision**, with the evidence
recorded, before any implementation.

## The decision

Input is **GFF3 + `annotation_summary.json`, not GenBank** — a candidate
that only reads GenBank is disqualified unless it also reads GFF.
Candidates from spec §8 item 4:

| Candidate | Note |
|---|---|
| pyCirclize (https://github.com/moshi4/pyCirclize) | Python, actively maintained, parses GFF as well as GenBank |
| plotgenes (https://github.com/aaronphillips7493/plotgenes) | From CLAW's author |
| Circos + custom GenBank→config converter | Heaviest; needs a format conversion we'd rather not own |
| **In-house Jinja + SVG** | Reads GFF directly; fits the existing report renderer container |

**Leaning in-house** (record the reasoning either way):
spec §6a.1 wants
clickable/hoverable features with tooltips carrying name, product,
strand, coordinates, feature source and task 40's
`pident`/`qcovhsp`/`bitscore`. A circos-style library emits a picture, not
an interactive DOM, so the tooltip layer would have to be rebuilt on top
of it anyway — and the in-house route adds no fourth container image,
reusing the C6/C7 report image
(spec §6a.5).
Weigh that against the effort of hand-rolling circular-genome layout.
Whichever wins, replace the `python:3.12-slim` + `TODO P4` placeholder in
`conf/containers.config`.

## Scope sketch

- Renderer reading GFF3 + `annotation_summary.json`; inline SVG output so
  it drops straight into task 43b's report.
- Per-feature tooltip metadata, including feature `source` (`miniprot` vs
  `mitos`) so a reader can separate the two without opening the summary
  JSON.
- **Dual-isoform rendering** (todo.md carry-forward): when
  `plastid_isoforms/` is present, walk **both `path1` and `path2`**, per
  spec §3.6 step 5.
  `ORGANELLE_MAP` is the *sole* stage aware of the plastid quadripartite
  structure — every other downstream stage consumes `target.fasta`
  unchanged.

## Sequencing note

Third of the four P4a tasks, deliberately. This is report *content*, not
a rendering blocker: task 43b can be built against the stub SVG and swap
in the real map when it lands. Do not let the tool decision block the
report tier.

## Outcomes

**Decision: in-house Jinja+SVG**, per the leaning above. All three
off-the-shelf candidates (pyCirclize, plotgenes, Circos) emit a static
picture with no DOM structure, so spec §6a.1's clickable/hoverable
per-feature tooltips would need a JS tooltip layer built on top of any
of them regardless — at which point the rendering step is the only
thing bought by adopting a dependency, and a circular-genome layout is
not a hard problem to hand-roll. The in-house route also reuses the
existing report/C6/C7 container instead of adding a fourth pinned image
(rule 12/19). `conf/containers.config`'s `ORGANELLE_MAP` entry now
points at the already-pinned `neoformit/daff-wf5-scripts` image (it
already carries `jinja2`, though the implementation ended up not
needing templating — a hand-built SVG string was simpler than a
template for this shape of output). Tooltips are a native SVG
`<title>` per feature — no JavaScript, matching the report's
self-contained-file requirement and the same choice task 47 makes for
the assembly graph. Spec §8 item 4 updated to record this as resolved.

**Two design gaps surfaced during implementation, not covered by the
original brief, and confirmed with the user before coding (both taken
as recommended):**

1. **Feature redundancy.** The merged annotation GFF is not
   deduplicated per locus — the broad miniprot pass (spec §8 item 3)
   hits every reference-panel ortholog at a locus independently, so a
   single physical gene routinely surfaces as many overlapping raw
   `mRNA`/`CDS` records (observed: 37 records for 13 true
   protein-coding loci on `INT-ANIMAL-01`). Drawing one arc per raw
   record would make the map an unreadable stack of near-duplicates.
   `bin/organelle_map.py`'s `cluster_loci()` collapses genomically-
   overlapping same-(seqid, strand) CDS calls onto one arc, keeping the
   best-scored (or, absent a score, longest-spanning) call as the
   primary annotation and recording every other gene label seen at
   that locus as `alt_genes` in the tooltip, so no candidate is
   silently dropped (CONSTITUTION principle 7). tRNA/rRNA calls have a
   single source (MITOS2) and are never clustered.
2. **Dual-isoform (`path2`) coordinate remap.** Spec §3.6 step 5
   requires `ORGANELLE_MAP` to render both plastid isoforms, but the
   annotation GFF only carries coordinates for `path1` (`=
   target.fasta`); `path2` differs only in the SSC segment being
   reverse-complemented. Implementing this required extending
   `plastid_canonicalise.py`'s `Result` (and hence
   `bin_metadata.json`'s `plastid_canonicalisation` block) with the
   individual `lsc_len`/`ir_len`/`ssc_len` segment lengths — previously
   only the edge *names* and total path lengths were recorded — and
   adding `bin_metadata.json` as a new `ORGANELLE_MAP` input (joined
   from `BIN_TARGET.out.metadata` in `main.nf`), which also solved a
   second latent gap: `bin_metadata.json` is the only place the map's
   true genome length lives (`target.fasta`'s length differs from any
   single contig's `length_bp` once the plastid substitution applies).
   `remap_to_path2()` mirrors any feature entirely inside the SSC
   segment and flips its strand; a feature spanning the LSC/SSC or
   SSC/IR boundary has no valid position in `path2`'s frame and is
   dropped from that panel only (with the primary `path1` panel
   unaffected).

**Scope note on "product".** The brief's tooltip field list included
`product`, but no stage in this pipeline stores a gene product
description — only short symbols (`COX1`, `trnG(tcc)`) are available
anywhere in the data model (`assets/loci.json` documents symbols only).
Fabricating product text was rejected as dishonest; the tooltip
presents the gene symbol as both name and product rather than
inventing a description field with nothing behind it.

**Files touched beyond the module itself:**
`bin/organelle_map.py` (new, 100% branch coverage via
`scripts/tests/test_organelle_map.py`), `bin/plastid_canonicalise.py`
(`Result` gains `lsc_len`/`ir_len`/`ssc_len`), `main.nf` (ORGANELLE_MAP
now joins `BIN_TARGET.out.metadata`), `conf/containers.config`, spec
§3.6 and §8 item 4.
