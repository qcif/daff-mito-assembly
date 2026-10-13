# Task 66 — Plastid map: quadripartite segment ring and exon structure

**Phase:** P4 (report content). Broken out of task
32_research_plastid_noncds.md §8, which evaluated the options and
recommended extending the in-house map (§8.3 option V1).

**Status:** open. **Best done after** task 63_cds_locus_clustering.md
(shared exon-overlap predicate, `rpl2` fix), task
64_protein_panel_label_check.md (removes the phantom `ycf3`/`cemA`
features) and task 65_plastid_noncds_barrnap.md (real tRNA/rRNA
features with plastid names). None is a hard precondition:
tRNA/rRNA features render from any GFF source.

## 1. Overview

The pipeline draws a circular genome map for every assembled organelle
(`ORGANELLE_MAP`, `bin/organelle_map.py`, task 44). Each annotated gene
is an arc on the circle, outer ring for the + strand and inner ring for
the − strand. Barcode genes get a printed label. Every feature gets a
hover tooltip naming it.

For a plastid (chloroplast) genome this map is missing the two things
a plant biologist reads first:

1. **The quadripartite structure.** A typical land-plant plastome is
   one circle made of four segments: a Large Single Copy region (LSC),
   a Small Single Copy region (SSC), and two identical Inverted Repeats
   (IRa, IRb) between them. Every published plastome map shows these
   segments as a ring, because they explain why genes like `rrn16` and
   `ndhB` appear twice. The pipeline already finds these segments while
   canonicalising the plastid assembly (spec §3.6) and records their
   lengths. The map just doesn't draw them.
2. **Exon structure.** About a dozen plastid protein-coding genes
   contain introns (`clpP`, `ycf3` and `rps12` have two each). The map currently draws
   each gene as one solid arc from its first base to its last, so
   intron sequence looks like coding sequence. The annotation already
   records each exon separately, so the map can show the exons as
   blocks joined by a thin line.

Both changes use data the stage already receives. Neither adds a tool,
container or pipeline input.

Gene labelling is **not** in scope. Non-barcode genes stay
tooltip-only, as on `animal_mt`; this was decided in task 32 §8.1.

## 2. Evidence from the fixture

From `INT-PLANT-01-pt` (task 32 §7.2, §8.1):

- `bin_metadata.json`'s `plastid_canonicalisation` block gives
  `branch: canonical`, `substitution_applied: true`, `lsc_len` 85,993,
  `ir_len` 25,490, `ssc_len` 18,303, `path1_len` 155,276. In path1
  coordinates that puts LSC at 1–85,993, IRb at 85,994–111,483, SSC at
  111,484–129,786 and IRa at 129,787–155,276. The segment layout is
  the same one `remap_to_path2()` already assumes.
- miniprot's exon calls agree with the `NC_018117.1` reference for
  `clpP` (3 exons), `ycf3` (3), `atpF`, `rpoC1`, `ndhB`, `rpl2` and
  `rps16`. These exons are real and worth drawing.
- miniprot can also produce absurd spliced alignments. Raw `cds.gff`
  contains one `rpl2` hit (`MP000542`) spanning 85,769–131,376, a
  45 kb "intron". `cluster_loci()` currently groups features by
  **span** overlap, so a call like this would absorb every
  same-strand gene inside its intron into one cluster. Today that is
  masked by the map drawing spans anyway. Once exons are drawn, the
  clustering basis matters.

## 3. Scope

### 3.1 Segment ring (`plant_pt` only)

- Draw a ring of four arcs **inside** the strand rings, one per
  segment. IRa and IRb share a colour; LSC and SSC each get their own.
  Use muted fills so the ring doesn't compete with gene features.
- Print the four segment names (`LSC`, `IRb`, `SSC`, `IRa`) on or next
  to their arcs. Four labels don't need collision handling.
- Each segment is its own SVG group with a `<title>` tooltip giving
  name, coordinates and length. Example tooltip:

  ```
  IRb (inverted repeat)
  path1:85994-111483
  25,490 bp
  ```

- Mark segment groups with their own attribute (e.g.
  `data-segment="IRb"`), **not** `data-gene`. Gene-oriented consumers
  select `g[data-gene]` and must not start treating segments as genes.
- **path2 panel:** same boundaries. The SSC tooltip and label should
  say it is inverted relative to path1, since that inversion is the
  only difference between the two isoforms.
- **Draw the ring only when the data supports it:**
  `plastid_path2_segments()` returns segment lengths (it already
  requires `substitution_applied` and three integer lengths), *and*
  `lsc + 2·ir + ssc` equals the drawn genome length. The same
  condition already triggers the path2 panel, so the ring and the
  second panel appear together. Consider renaming the helper
  (e.g. `plastid_segments()`) now that it serves both. Otherwise the output is
  unchanged. That covers the non-canonical branch, the withheld
  substitution, `plant_mt` and `animal_mt`. Never draw a guessed
  segment ring.
- Add the segment colours to the legend, but only on maps that draw
  the ring.

### 3.2 Exon blocks (all targets)

- Carry the per-exon coordinates that `parse_annotation_gff()` already
  returns (`cds` list) through `cds_features()` into each locus.
- Render a multi-exon locus as one group (`g[data-gene]`, one tooltip)
  containing one wedge per exon plus a thin arc connector along the
  middle of the strand ring between consecutive exons.
- Add the exon count and exon coordinates to the tooltip when there is
  more than one exon:

  ```
  clpP
  CDS | path1:71901-73900 (-) | 3 exons
  exons: 71901-72126, 72750-73042, 73830-73900
  source: miniprot
  ...
  ```

- A single-exon locus must render **byte-identically** to today. This
  keeps the change a no-op on almost all of `animal_mt`.
- tRNA/rRNA rows are single features in every annotator we use
  (MITOS2, barrnap); leave them as single wedges.
- `remap_to_path2()` must mirror every exon of an SSC-internal locus,
  not just the outer span. Exon order should follow increasing
  coordinate in the new frame. The existing rule still applies: a
  locus with any exon straddling an SSC boundary is dropped from the
  path2 panel only.

### 3.3 Cluster on exon overlap, not span

- Change `cluster_loci()`'s overlap test so two calls cluster only
  when their **exon blocks** overlap on the same seqid and strand.
  Span overlap alone (one call sitting in the other's intron) is not
  enough. Import the exon-overlap predicate task
  63_cds_locus_clustering.md §2.1 adds to a shared `bin/` module.
  Don't write a third copy.
- Primary selection inside a cluster (`_score_key`: best score, else
  longest) is unchanged.
- Keep "longest span" as the fallback tiebreak, but check that a
  pathological long-intron call (the 45 kb `MP000542` shape) cannot
  win on span alone over a well-scored compact call at the same
  locus. If it can, switch the fallback to total exon length and
  record it in Outcomes.

### 3.4 Report

- `scripts/report/templates/components/assembly.html` promotes native
  `<title>` tooltips to Bootstrap tooltips via
  `promoteSvgTitleTooltips('organelle-map', 'g[data-gene]', ...)`.
  Extend the selector, or add a second call, so segment tooltips get
  the same treatment in the inline map and in the fullscreen modal.
- The standalone `<id>.map.svg` keeps native `<title>` tooltips
  (spec §6a), so it needs no change beyond §3.1–3.3.

### 3.5 Out of scope

- Printed labels for non-barcode genes (decided: tooltip only).
- Zoom, and any JS mapping library (CGView.js: task 32 §8.2–8.3
  option V2, deferred).
- GC-content ring and functional-group colouring (task 32 §8.1
  nice-to-haves; would need `target.fasta` as a new input).
- Normalising barrnap names (`tRNA-Met` → `trnM-CAU`,
  `16S_rRNA` → `rrn16`). That is done at ANNOTATE merge time by task
  65_plastid_noncds_barrnap.md.
- `rps12` *trans*-splicing. miniprot never calls the 5′ exon (task 32
  §7.2), so there is nothing to join.
- The ANNOTATE merge dropping `rpl2`'s second IR copy. That is a
  merge bug, fixed in task 63_cds_locus_clustering.md.

## 4. Pseudocode

```
segments = plastid_path2_segments(metadata)  # (lsc, ir, ssc) or None
if segments and lsc + 2*ir + ssc == length:
    ring = [("LSC", 1, lsc), ("IRb", lsc+1, lsc+ir),
            ("SSC", lsc+ir+1, lsc+ir+ssc), ("IRa", lsc+ir+ssc+1, length)]
    draw ring at radius < inner strand ring; label + tooltip each arc

for locus in clustered_loci:                  # clustering on exon overlap
    exons = locus.exons or [(locus.start, locus.end)]
    if len(exons) == 1: draw exactly as today
    else:
        open <g data-gene=...><title>…exons…</title>
        for (a, b) in exons: wedge(a, b)
        for gap in consecutive pairs: thin connector arc at ring mid-radius
        close </g>
```

## 5. Testing

### 5.1 Unit tests (`scripts/tests/test_organelle_map.py`, CONSTITUTION rule 14/15)

Keep `bin/organelle_map.py` at its current 100 % branch coverage. Add
cases for:

- **Segment ring present:** canonical metadata with
  `substitution_applied: true` gives four `data-segment` groups with
  the expected coordinates. IRa and IRb share a fill. The legend gains
  the segment entries.
- **Segment ring absent:** withheld substitution, non-canonical
  branch, missing lengths, non-integer lengths, lengths that don't sum
  to the genome length, and a non-plastid target. Each gives no
  `data-segment` and no legend change.
- **path2 ring:** same boundaries, and the SSC tooltip notes the
  inversion.
- **Single-exon locus is byte-identical** to the pre-change output.
  Pin this with a fixed expected string or a hash.
- **Multi-exon locus:** one `g[data-gene]`, N wedges, N−1 connectors,
  and a tooltip listing the exons.
- **path2 remap of a multi-exon SSC locus:** exons mirrored, reordered
  and strand flipped. A locus with one exon straddling the SSC boundary
  is dropped from path2 and kept on path1.
- **Clustering:** a call nested inside another's intron stays a
  separate locus. Calls with overlapping exons still cluster. A long
  span with poor or no score does not displace a compact well-scored
  call (§3.3).

Mock nothing below the module boundary; the inputs are small GFF and
JSON files written in `tmp_path` (task 27's boundary-mocking rule).

### 5.2 Stub run

No channel or process signature changes, so `-stub-run` should pass
unchanged. Confirm it.

### 5.3 Integration

`tests/integration/assertions.sh` currently checks the assembly-graph
SVG (task 47) but not the organelle map. Add a light structural check
in the same style:

- for every assembling sample, `annotation/<id>.map.svg` exists,
  contains `<svg` and at least one `data-gene=`;
- for `INT-PLANT-01-pt`, it also contains four `data-segment=` groups
  per panel;
- for `INT-ANIMAL-01`, it contains no `data-segment=`.

No new fixtures are needed. `INT-PLANT-01-pt` exercises the canonical
branch and has multi-exon genes; `INT-ANIMAL-01` exercises the
no-ring, mostly single-exon path.

### 5.4 Visual check

Render both fixtures' reports and inspect the map by eye, inline and
in the fullscreen modal. Give the report paths to the user, as
CLAUDE.md requires for notable report changes.

## 6. Spec updates

- spec/06a-reports.md, Assembly tab, map description: add the segment
  ring on `plant_pt` (drawn only when canonicalisation substituted
  path1) and exon-block rendering. State that non-barcode genes are
  tooltip-only on every target.
- spec/03-organelles.md §3.6 step 5: note that `ORGANELLE_MAP` also
  draws the LSC/IR/SSC ring from the recorded segment lengths, not only
  the path2 remap.
- spec/07-open-questions.md §8 item 4 (organelle map renderer, resolved):
  one line noting the plastid extension, and that CGView.js was
  evaluated and deferred in task 32_research_plastid_noncds.md.

## 7. Checklist

- [ ] §3.1 segment ring, with all degrade paths
- [ ] §3.2 exon blocks, single-exon byte-identical, path2 remap
- [ ] §3.3 exon-overlap clustering and the long-intron tiebreak check
- [ ] §3.4 report tooltip binding for segments
- [ ] §5.1 unit tests, 100 % branch coverage held; flake8 clean
- [ ] §5.2 stub run passes
- [ ] §5.3 integration assertions added and passing
- [ ] §5.4 reports rendered and paths given to the user
- [ ] §6 spec updated
- [ ] Outcomes section written
