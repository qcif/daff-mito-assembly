# Task 44 — Stage 14 `ORGANELLE_MAP`: annotated organelle diagram

> **STUB — expand before executing.**

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
