# Task 19 — Plastid canonicalisation algorithm spec (C4 clean-room read)

**Phase:** P3 (from [spec §6](../../spec/06-phases.md)).
**Goal:** Read
[`reference-material/ptgaul/combine_gfa.py`](../../reference-material/ptgaul/combine_gfa.py)
and — combined with existing spec material
([spec §3.6](../../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation))
— produce `spec/plastid-canonicalisation.md`: a self-sufficient
algorithm design document that will serve as the sole permitted
reference for the C4 implementation task
([task 20](20_plastid_canonicalise.md)).

**No runnable code is produced by this task.** The deliverable is a
prose + pseudocode + test-case-matrix document. Successful completion
means an implementer can, without ever opening `combine_gfa.py`, write
`bin/plastid_canonicalise.py` from `spec/plastid-canonicalisation.md`
alone.

---

## Why this task exists — clean-room firewall

[`reference-material/ptgaul/combine_gfa.py`](../../reference-material/ptgaul/combine_gfa.py)
ships **no licence**; the ptGAUL README declares no terms. Default
copyright ("all rights reserved") applies. Copying the source into
this repo is prohibited.

The algorithm itself — a bioinformatics recipe: longest edge → LSC,
deepest edge → IR, remainder → SSC; concatenate to form two isoforms —
is not copyrightable. But **implementing an algorithm without
inadvertently copying its expression** is easier to defend when the
person reading the source and the person writing the replacement are
separated by a written specification.

That is what this task produces: the specification that acts as the
firewall between "read ptGAUL" (this task) and "write Python"
([task 20](20_plastid_canonicalise.md)).

**Constraint on the executor of this task:**

- **You may open `combine_gfa.py`.** In fact you must — you are the
  reader.
- **You must not write any runnable code as a deliverable of this
  task.** Pseudocode in the spec is fine; committing
  `bin/plastid_canonicalise.py` is out of scope.
- Anything that reads as "please copy this identifier verbatim into
  Python" is out of scope — the spec describes *what* the code must
  compute, not the exact identifiers, argument names, or comment
  wording ptGAUL uses.

**Constraint on the executor of task 20 (documented here for
context):** task 20 must not open `combine_gfa.py` and must not read
this task's Outcomes section if it references ptGAUL constructs
directly. The finalised `spec/plastid-canonicalisation.md` is task
20's only permitted algorithm reference.

---

**Prerequisite:** None. This task is a documentation-only exercise
against static reference material that already exists in the repo.

**Exit criteria:**

- [`spec/plastid-canonicalisation.md`](../../spec/plastid-canonicalisation.md)
  exists, committed to the repo, and covers the sections listed in
  §2 below.
- The document is *self-sufficient* for implementation: no `TODO`,
  no `see combine_gfa.py`, no `check ptGAUL for …` references. Every
  detail the implementer needs is either resolved in-doc or
  explicitly flagged as an implementation-choice-left-to-C4 with
  guidance on the acceptable range of choices.
- No new files created under `bin/`, `scripts/tests/`, or
  `modules/local/`. No changes to `nextflow.config`,
  `conf/containers.config`, or `main.nf`.
- The document renders cleanly (markdown lint if configured; visual
  inspection otherwise).
- Task 18 and downstream cross-references that point at
  `spec/plastid-canonicalisation.md` resolve
  (e.g., [task 20 §1](20_plastid_canonicalise.md#1-binplastid_canonicalisepy-new-custom-logic-c4)
  links here).

**Not in scope:**

- Any implementation work — `bin/plastid_canonicalise.py`, C3 hook,
  isoforms emit, unit tests, integration assertions, bin_bounds
  update. All of that lives in [task 20](20_plastid_canonicalise.md).
- Container work. Uses no containers.
- Modifying spec §3.6. It stays as the pipeline-flow summary; the
  new doc is the deep-dive companion (following the same
  spec-file-plus-adjunct-doc pattern used by
  [`spec/miniprot.md`](../../spec/miniprot.md)).
- Deciding whether the ptGAUL algorithm is the "right" plastid
  canonicalisation approach. That decision is already made
  ([spec §3.6](../../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation));
  this task documents *how* to implement it.

---

## 1. Reading protocol

1. Read `combine_gfa.py` end-to-end. Take notes in a scratch buffer
   *outside* the repo (do not commit working notes — see §4 below).
2. Read [spec §3.6](../../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation)
   in the current spec to understand the algorithm's role in the
   pipeline and the existing algorithmic sketch.
3. Read [`reference-material/ptgaul/ptGAUL.sh`](../../reference-material/ptgaul/ptGAUL.sh)
   only if needed to understand how `combine_gfa.py`'s inputs are
   produced upstream (edge FASTA, sorted depth TSV) — this matters
   because our C4 will parse GFA directly rather than accept
   pre-processed inputs.
4. Read Flye's own GFA output format
   ([Flye user guide](../../reference-material/flye-user-guide.md))
   to nail down the `S`-line depth tag convention (`dp:f:`).

## 2. Deliverable — `spec/plastid-canonicalisation.md` structure

Author the document to the following outline. Length target: ~200–400
lines. Longer is fine if edge cases warrant it; shorter probably
means an implementer will end up asking questions.

1. **Purpose and pipeline context** — one paragraph, cross-links to
   [spec §3.6](../../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation)
   and [spec §2 stage 10 (BIN_TARGET)](../../spec/02-stages.md#2-stage-detail).
2. **Biological background** — quadripartite plastid structure (LSC /
   IR / SSC), why the IR appears at 2× coverage, why the two SSC
   orientations are both valid. Two paragraphs max — this exists to
   ground the algorithm, not to teach botany.
3. **Input specification** —
   - GFA file layout (which lines matter, which are ignored).
   - `S` line schema, `dp:f:` / `DP:f:` depth tag, missing-tag
     behaviour, sequence-in-column-3 assumption.
   - Explicitly document what `L` (link) lines are and why the
     algorithm ignores them.
   - Assumptions the caller must guarantee (well-formed Flye GFA,
     single connected component or not — pick a stance and justify).
4. **Output specification** —
   - `path1.fasta`, `path2.fasta` — record IDs, header shape, single
     record per file, uncompressed FASTA.
   - Metadata dict shape: field names, types, permitted values for
     `branch`, what null / not-applicable looks like on non-canonical
     branches.
   - Where they land in the process work-dir (`outdir/`); how the
     module publishes the `plastid_isoforms/` directory.
5. **Branch classification rules** — decision table:

   | Edge count | Branch | Output | Downstream effect |
   |---|---|---|---|
   | 1 | `resolved_circle` | none | C3 keeps its single-edge target unchanged |
   | 3 | `canonical` (provisional) | proceed to §6 | If §6 does not degenerate, substitute path1 into target.fasta |
   | anything else | `non_canonical` | none | diagnostic marker in bin_metadata.json |

   Also enumerate the degenerate sub-cases inside the 3-edge branch
   (e.g., longest and deepest are the same edge, missing depth
   tags → all zero, etc.) and specify the fall-back to
   `non_canonical`.
6. **Canonicalisation algorithm (3-edge)** — pseudocode showing:
   - Selection: longest → LSC; deepest → IR; the remaining edge → SSC.
     Explicit tie-breaker rules.
   - Isoform construction: `path1 = LSC + IR + SSC + rc(IR)`;
     `path2 = LSC + IR + rc(SSC) + rc(IR)`.
   - Sanity-check invariant: `len(path1) == len(path2) == LSC + 2·IR + SSC`.
   - Why `path1` uses SSC and `path2` uses `rc(SSC)` (biological
     reason, one sentence).
7. **CLI arg-surface requirements** — describe the shape task 20
   should design (positional GFA, `--outdir`, `--json-out`), and
   explicitly say *not* to mirror `combine_gfa.py`'s `-e / -d / -o`
   triple. State that the CLI is for debugging only; the library
   entry point is the contract.
8. **Library entry point requirements** — signature, return type
   shape (e.g., a `Result` namedtuple with the fields listed in
   §4 metadata), side effects (writes files iff `canonical`).
9. **Test-case matrix** — a table of at least the following synthetic
   GFAs and expected outcomes:
   - 3-edge canonical (LSC 90k depth 1, IR 25k depth 2, SSC 15k
     depth 1) → `canonical`.
   - 3-edge with `DP:f:` (upper-case) depth tag → still parses.
   - 3-edge with missing depth tag on the true IR → falls to
     `non_canonical` (safety).
   - 3-edge with LSC and IR colliding on the same edge → falls to
     `non_canonical` (safety).
   - 1-edge resolved circle → `resolved_circle`.
   - 2-edge → `non_canonical`.
   - 4-edge → `non_canonical`.
   - Round-trip length invariant on the canonical case.
   - `path2` uses `rc(SSC)` at the correct offset.

   Each row should include the *inputs* the test constructs and the
   *assertions* it makes — enough for task 20 to write the tests
   without design decisions.
10. **Implementation-choice boundaries** — anywhere the algorithm
    leaves latitude (e.g., which biopython class to hold the
    sequence, whether the depth regex is greedy, how empty
    directories are handled), state the acceptable range explicitly.
    The implementer should never need to guess.
11. **Provenance** — one paragraph acknowledging ptGAUL as the
    algorithm's original expression and this repo's clean-room
    posture. Cite:
    - [`reference-material/ptgaul/combine_gfa.py`](../../reference-material/ptgaul/combine_gfa.py)
    - the ptGAUL paper (if the executor can find a citation) or the
      GitHub URL as fallback
    - [CONSTITUTION.md rule 12](../../CONSTITUTION.md) (biocontainer
      preference; explains why a rewrite is preferred over vendoring)

## 3. Cross-references to add

Once `spec/plastid-canonicalisation.md` is committed:

- Add a link from
  [spec §3.6](../../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation)
  ("See [spec/plastid-canonicalisation.md](../../spec/plastid-canonicalisation.md)
  for the full algorithm specification.") — one line, no other §3.6
  changes.
- Confirm [task 20 §1](20_plastid_canonicalise.md#1-binplastid_canonicalisepy-new-custom-logic-c4)
  link to the new doc resolves.

## 4. Working-notes hygiene

- Take reading notes outside the repo. Do not commit a
  `notes/` directory, a scratch file under `tmp/`, or annotations
  inside `reference-material/`.
- If ptGAUL identifiers or code fragments appear in your notes, do
  not carry them into `spec/plastid-canonicalisation.md`. The spec
  should read as if the author has never seen `combine_gfa.py` — it
  describes the algorithm in the language of the pipeline
  (`assembly_target`, `plant_pt`, `bin_metadata`), not ptGAUL's
  (`myList`, `SSR_seq`, `sorted_depth_file`).
- The Outcomes section (§7 below) may cite ptGAUL constructs *as
  provenance metadata* but not as implementation guidance. Task 20
  is instructed to treat Outcomes as provenance-only.

## 5. Verification

1. `spec/plastid-canonicalisation.md` exists and covers all sections
   listed in §2.
2. `grep -n "TODO\|FIXME\|see combine_gfa\|check ptGAUL" spec/plastid-canonicalisation.md`
   returns no matches.
3. [`spec/03-organelles.md`](../../spec/03-organelles.md) §3.6 links to
   the new doc.
4. [task 20](20_plastid_canonicalise.md)'s links to
   `spec/plastid-canonicalisation.md` resolve.
5. Sanity read: a fresh implementer, given only the new spec and
   this repo's normal onboarding docs, could plausibly implement
   `bin/plastid_canonicalise.py`. If not, the spec is under-detailed
   and needs more work.

## 6. Deliverables checklist

- [x] [`spec/plastid-canonicalisation.md`](../../spec/plastid-canonicalisation.md)
      created with all sections from §2 above.
- [x] [`spec/03-organelles.md`](../../spec/03-organelles.md) §3.6 has a
      one-line link to the new doc (was already present from prior spec
      authoring; no change needed).
- [x] No new code files under `bin/`, `scripts/tests/`, or
      `modules/local/`.
- [x] No changes to `nextflow.config`, `conf/*.config`, or `main.nf`.
- [x] Working notes are not committed.
- [x] `grep -n "TODO\|FIXME" spec/plastid-canonicalisation.md`
      returns nothing.

## 7. Notes / non-issues

- **This is not a "write a design doc" boilerplate exercise.** The
  quality of the spec directly determines the risk profile of task
  20. Under-specify and task 20's implementer will either (a) file
  a follow-up back to task 19 (fine, slow) or (b) reach for
  `combine_gfa.py` (bad, defeats the point). Aim for (c): the spec
  is a joy to implement from.
- **Ambiguity is fine, ambiguity-in-hiding is not.** If a detail
  genuinely admits multiple valid choices (e.g., handling of a
  malformed `S` line), state the choices and let task 20 pick. Do
  not silently omit the decision.
- **Spec author == task 19 executor.** For this task the "reader"
  and "spec writer" are the same person; the firewall is between
  this task and task 20. In an ideal world these would be two
  people. Enforcement is procedural (task-file separation) rather
  than personnel-based.
- **If two spec versions are being iterated on, only the committed
  one is authoritative.** Task 20's implementer must not read
  draft/branch versions of this spec. That is why the deliverable
  is a single committed file, not a PR-embedded draft.

## 8. Outcomes

- Path to committed spec: `spec/plastid-canonicalisation.md`
  (commit `6a959ac`)
- Word count / line count: 455 lines
- ptGAUL constructs cited here for provenance (not to be used by
  task 20):
  - The ptGAUL depth-sorting step used
    `awk -F ":" '/^S/{print substr($1,3,7)"\t"$3}'` on the GFA,
    relying on `dp:f:` being parseable as the third `:` -delimited
    field. Our spec mandates a more robust regex scan instead.
  - ptGAUL variable names (`myList`, `SSR_seq`, `sorted_depth_file`,
    `longest`, `IR`) appear nowhere in the delivered spec.
  - `combine_gfa.py` does not guard against degenerate cases (all-equal
    depth, LSC/IR collision, missing tags). The spec adds these
    degenerate checks as explicit `non_canonical` branches.
- Ambiguities left open for the implementer with rationale:
  - **Tie-breaking rule**: when two edges share the maximum length or
    depth, the spec requires *a* deterministic rule but does not
    mandate lexicographic order; any documented deterministic rule is
    acceptable. Rationale: this choice has no biological significance
    and the implementer can pick what reads most clearly.
  - **`LN:i:` vs `len(sequence)` for edge length**: spec mandates
    `len(sequence)`. Rationale: LN tag may be absent in minimal
    synthetic GFAs used in tests; the sequence is always present.
  - **Malformed depth float** (e.g., `dp:f:nan`): either treat as 0.0
    or raise ValueError. Rationale: this scenario is unlikely from Flye
    and either choice is safe.

## 9. Provenance and licensing (relocated from the spec)

This section is the project's **single record** of where the
canonicalisation algorithm originated. It was moved here out of
`spec/plastid-canonicalisation.md` so that the specification — the
task 20 implementer's only permitted reference — carries no pointer to
upstream source. Do not re-introduce this material into `spec/`.

- Original expression of the algorithm:
  `reference-material/ptgaul/combine_gfa.py`, part of the ptGAUL
  pipeline (GitHub: `https://github.com/Bean061/ptGAUL`). Retained
  under `reference-material/` for provenance only — never copied,
  imported, or vendored.
- Paper: Xu L, Dong Z, Fang L, et al. (2023). "ptGAUL: A pipeline for
  the assembly and classification of plant organellar genomes."
  *Molecular Ecology Resources*. Cite the GitHub URL as fallback.
- **Licensing finding:** the ptGAUL repository ships no licence file
  and its README declares no terms, so default copyright ("all rights
  reserved") applies to its source code. The algorithm itself — a
  recipe for identifying quadripartite regions by sequence length and
  read depth — is not copyrightable expression, so re-implementing it
  from a prose specification is sound while copying the source is not.
- **Consequence:** `bin/plastid_canonicalise.py` is a clean-room
  re-implementation. This task read the upstream source and wrote the
  specification; task 20 reads only the specification and writes the
  code. Neither task touches both sides of that boundary.
  Note: CONSTITUTION.md has **no** rule covering vendoring or
  licensing — an earlier draft of the spec cited "rule 12" for this,
  which is actually about biocontainers. If this policy should be
  binding project-wide it needs to be added to the constitution
  properly; tracked in [todo.md](../todo.md).
