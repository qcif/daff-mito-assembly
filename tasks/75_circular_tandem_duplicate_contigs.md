# Task 75 — Circular contigs that are a tandem-duplicated genome (stub)

**Status:** stub. Scope to be decided when picked up.

**Related spec:** spec/03-organelles.md §3.1 (typical sizes per target),
§3.7 (BIN_TARGET), §3.7.4 (circularity comes from Flye). CONSTITUTION
principle 7 (negative clarity).

## Overview

metaFlye can close a circular organelle contig that holds the genome
**twice** in a row (a tandem duplicate, or concatemer) instead of once.
The contig is flagged circular, aligns well to the organelle panel and
carries every gene, so BIN_TARGET selects it and the barcodes extracted
from it look normal. But its length is about double the real genome,
and anything that relies on the assembled sequence as a faithful
genome (length, gene order, the organelle map, the deliverable FASTA)
is wrong.

Found in task 68_animal_mt_microsatellite_recruitment.md (2026-10-20).
On the `INT-ANIMAL-01` fixture (*Acyrthosiphon pisum*, mitogenome
about 17 kb):

- Recruited against `v2026.09_1`: the target is one 16,799 bp circular
  contig.
- Recruited against `v2026.10` (repeat-masked `animal_mt` panel):
  the target is one 33,905 bp circular contig (Flye `circ.=Y`,
  coverage 58x, multiplicity 9). 82 % of the 31-mers in its first
  17 kb also appear in the rest of the contig.

The only change is a smaller recruited read pool (5.57 Mb to 4.99 Mb),
which changes Flye's subsample above `max_cov`. So this is a property of
the assembler on a given read set, not of the masking itself, and it
can happen on real samples too. Today it is only caught by the
integration test's target-length bound; on a client sample nothing
would flag it.

## Questions to answer when scoped

- How often does it happen on real data? Check the target contig
  length against the expected size range (spec §3.1) across past
  remote and Azure runs.
- Detection: a self-alignment or k-mer test for a contig that repeats
  itself over most of its length. Where does it belong (BIN_TARGET,
  after selection)?
- Response: trim to one copy (and re-check circularity and gene
  order), or flag only (a warning and a report finding) and leave the
  sequence untouched? Trimming changes the deliverable and needs
  explicit sign-off.
- Does Flye have a setting that avoids this (for example the
  `--meta` / subsampling interaction), and is the subsample seed
  pinned? See tasks/todo.md, Recruitment: `INT-ANIMAL-01`'s assembly is
  not stable across recruitment changes.
- Integration fixture: a pinned read set that reproduces the doubled
  contig would make a good regression test.

The `INT-ANIMAL-01` fixture was rebuilt in task 68 so it no longer
reproduces the doubled contig; a regression fixture for this task would
need to be built deliberately.

Python changes need unit tests in `scripts/tests/` (CONSTITUTION
rule 14).
