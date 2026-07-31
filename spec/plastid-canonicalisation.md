# Plastid Canonicalisation Algorithm Specification

## 1. Purpose and pipeline context

This document is the authoritative algorithm specification for
`bin/plastid_canonicalise.py` (custom-logic component C4 in
[spec §2.2](02-stages.md#22-custom-logic-components)). It is produced
by task 19 and is the **sole permitted algorithm reference** for the
task 20 implementation. After reading this document an implementer must
be able to write `bin/plastid_canonicalise.py` without consulting any
other source, including `reference-material/ptgaul/combine_gfa.py`.

C4 is invoked from within C3 (`bin/bin_target.py`) on the `plant_pt`
branch of stage 10 (BIN_TARGET); see
[spec §3.6](03-organelles.md#36-plastid-quadripartite-canonicalisation-ptgaul-derived)
for the pipeline-flow summary and
[spec §2 stage 10](02-stages.md#2-stage-detail) for BIN_TARGET's full
role. C4 is not a new pipeline stage; it is a helper that runs inside
the existing BIN_TARGET process and shares its container.

---

## 2. Biological background

Land-plant plastid genomes have a canonical **quadripartite** structure:
a large single-copy region (LSC, ~80–90 kb), an inverted repeat A
(IRa, ~20–30 kb), a small single-copy region (SSC, ~15–20 kb), and an
inverted repeat B (IRb), where IRb is the reverse complement of IRa. The
two inverted repeats flank the SSC. Because the IR is present twice in
the plastome, sequencing reads derived from the IR accumulate at roughly
double the per-base coverage of the single-copy regions.

The SSC region can insert in either orientation relative to the LSC and
IR — both orientations are biologically real and co-exist within a cell.
A plastid population is therefore a 50/50 mixture of two isoforms. A
"complete" plastid assembly must represent both. The two isoforms differ
only at the SSC segment: one uses the SSC sequence as assembled, the
other uses its reverse complement.

---

## 3. Input specification

### 3.1 File format

The sole input is a Flye `assembly_graph.gfa` file. The GFA 1 format
uses tab-separated lines, each beginning with a single-character record
type. Only two record types are relevant:

- **`S` lines** (Segment): carry edge name, sequence, and optional
  tag fields including depth. These are the only lines the algorithm
  reads.
- **`L` lines** (Link): describe connections between segment ends.
  The canonicalisation algorithm ignores `L` lines entirely. Edge
  connectivity is irrelevant because the quadripartite structure is
  recovered solely from length and depth.
- All other record types (`H`, `P`, `W`, etc.) are ignored.

### 3.2 S-line schema

A Flye `S` line has the following tab-separated columns:

```
S    {edge_name}    {sequence}    [optional_tag …]
```

Columns:
1. `S` — record type character.
2. Edge name — a string such as `edge_000001`. This is the unique
   identifier for the segment within this GFA file.
3. Sequence — the DNA string (A/C/G/T, upper-case). Flye does not emit
   `*` (absent sequence) for organelle assemblies; the implementer may
   treat `*` as a fatal input error.
4. Optional tag fields — one or more `TAG:TYPE:VALUE` triplets separated
   by tabs. The only tag of interest is the depth tag.

### 3.3 Depth tag

Flye writes the per-edge mean read depth as:

```
dp:f:{value}
```

where `dp` is the tag name, `f` denotes a float, and `{value}` is a
non-negative floating-point number (e.g., `dp:f:56.7`).

**Case sensitivity:** Some GFA-producing tools emit `DP:f:` (upper-case
tag name). The parser must accept both `dp:f:` and `DP:f:`. A
case-insensitive regex match on the tag name is the recommended approach.

**Tag position:** The depth tag may appear at any position among the
optional fields. Do not assume it is the first or second optional field;
scan all optional fields for a match.

**Missing depth tag:** If an `S` line has no depth tag at all (neither
`dp:f:` nor `DP:f:`), treat that edge's depth as **0.0**. An edge with
depth 0.0 will never be selected as IR (deepest) unless all other edges
also have depth 0.0, in which case the assembly falls to `non_canonical`
via the degenerate-depth check in §5.3.

### 3.4 Sequence length

Compute edge length as `len(sequence)` from column 3, not from a
`LN:i:` tag. The sequence is always present; the tag may be absent or
disagree on truncated test inputs.

### 3.5 Caller guarantees and assumptions

The algorithm assumes:

- The GFA was emitted by Flye's metagenomic assembler (`flye --meta`)
  and has passed through the METAFLYE and BIN_TARGET filtering steps.
- All `S` line sequences are uppercase DNA strings.
- The GFA file is well-formed (parseable line by line; no binary
  content).
- The algorithm does **not** require the GFA to form a single connected
  component. Link topology is ignored entirely.
- Edge names within a single GFA file are unique. If duplicate edge
  names are encountered, raise a fatal `ValueError` — this indicates
  a malformed or concatenated input.

---

## 4. Output specification

### 4.1 FASTA files

On the `canonical` branch (3-edge case that passes degenerate checks),
C4 writes exactly two files into `{outdir}/`:

| File | Record ID | Content |
|---|---|---|
| `path1.fasta` | `path1` | LSC + IR + SSC + rc(IR) |
| `path2.fasta` | `path2` | LSC + IR + rc(SSC) + rc(IR) |

Each file contains **one FASTA record**: a header line (`>{id}`) followed
by the complete concatenated sequence on a single line (no line-wrapping).
Files are uncompressed plain text. The record IDs are the literal strings
`path1` and `path2`; they do not reference the original edge names.

No output files are written on the `resolved_circle` or `non_canonical`
branches. The `outdir/` directory must be created by the caller before
invoking C4, or by C4 itself (using `os.makedirs(outdir, exist_ok=True)`
— implementation choice).

### 4.2 Metadata dict (return value)

C4 returns a `Result` namedtuple (or equivalent typed dict) with the
following fields:

| Field | Type | Description |
|---|---|---|
| `branch` | `str` | One of: `"canonical"`, `"resolved_circle"`, `"non_canonical"` |
| `edge_count` | `int` | Number of `S` lines parsed from the GFA |
| `lsc_edge` | `str \| None` | Edge name assigned LSC; `None` on non-canonical branches |
| `ir_edge` | `str \| None` | Edge name assigned IR; `None` on non-canonical branches |
| `ssc_edge` | `str \| None` | Edge name assigned SSC; `None` on non-canonical branches |
| `path1_len` | `int \| None` | Length of path1 sequence in bases; `None` if no paths written |
| `path2_len` | `int \| None` | Length of path2 sequence in bases; `None` if no paths written |
| `non_canonical_reason` | `str \| None` | Short reason string on `non_canonical` branch (see §5.3); `None` otherwise |

The caller (`bin/bin_target.py`) merges this dict into `bin_metadata.json`
under the key `"plastid_canonicalisation"`.

### 4.3 Published isoform directory

The `BIN_TARGET` Nextflow module publishes `${outdir}/path1.fasta` and
`${outdir}/path2.fasta` as `plastid_isoforms/path1.fasta` and
`plastid_isoforms/path2.fasta` in the sample's output bundle. This
publish step is wired in the module definition (task 20 scope), not
inside C4.

---

## 5. Branch classification rules

### 5.1 Primary decision table

After parsing all `S` lines and counting edges:

| Edge count | Branch | C4 output | Downstream effect |
|---|---|---|---|
| `1` | `resolved_circle` | No files written | C3 keeps the single-edge sequence as `target.fasta` unchanged; no isoforms emitted |
| `3` | Proceed to §5.2 degenerate checks | Potentially `canonical` | If §5.2 passes: write path1/path2, substitute path1 into target.fasta |
| Any other count | `non_canonical` | No files written | Diagnostic marker written to bin_metadata.json; BANDAGE_NG output is the operator's review tool |

### 5.2 3-edge degenerate checks

Before applying the canonicalisation algorithm, the 3-edge case must
pass all of the following checks. Any failed check causes C4 to return
`non_canonical` immediately (before edge assignment), with an
appropriate `non_canonical_reason`:

1. **Depth collision check** — the edge selected as "deepest" (IR
   candidate) must be strictly deeper than at least one other edge.
   If all three edges share the same depth value (including the case
   where all depth tags are absent → all 0.0), classification is
   ambiguous; fall back to `non_canonical`.
   Reason string: `"depth_tie: all edges have equal depth"`.

2. **LSC–IR identity check** — the "longest" edge (LSC candidate) and
   the "deepest" edge (IR candidate) must not be the same edge. If
   length-maximum and depth-maximum are the same edge, the assembly is
   ambiguous.
   Reason string: `"lsc_ir_collision: longest and deepest are the same edge"`.

3. **Length degenerate check** — an edge of zero length (empty sequence
   string) is an invalid input. If any edge has length 0, return
   `non_canonical`.
   Reason string: `"zero_length_edge"`.

### 5.3 Non-canonical reason strings

| Situation | `non_canonical_reason` |
|---|---|
| Edge count ≠ 1 and ≠ 3 | `"edge_count:{n}"` where `{n}` is the actual count |
| All depths equal (including all-zero) | `"depth_tie: all edges have equal depth"` |
| LSC and IR would be the same edge | `"lsc_ir_collision: longest and deepest are the same edge"` |
| Any edge has zero-length sequence | `"zero_length_edge"` |

The `non_canonical_reason` field is `None` on `canonical` and
`resolved_circle` branches.

---

## 6. Canonicalisation algorithm (3-edge)

This section applies only when the 3-edge case has passed all
degenerate checks in §5.2.

### 6.1 Edge assignment

```
edges = list of (name, sequence, depth) triples parsed from all S lines

LSC_edge  = edge with maximum len(sequence)    # longest → LSC
IR_edge   = edge with maximum depth            # deepest → IR
SSC_edge  = the one remaining edge
```

**Tie-breaking rules:**

- **Length tie** (two edges share the maximum length): select the one
  with the higher depth as LSC. If depth is also tied, select the one
  whose edge name sorts first lexicographically. This is an
  implementation-choice boundary (§10); the spec requires a
  deterministic rule, not a specific one.
- **Depth tie at the top** (two edges share the maximum depth): select
  the one with the greater sequence length as IR. If length is also
  tied, select the one whose edge name sorts last lexicographically.
  Again, determinism is the requirement.
- After the depth-tie, the LSC–IR identity check (§5.2 item 2) may
  still fire; in that case the result is `non_canonical` as documented.

### 6.2 Isoform construction

```
IR_rc  = reverse_complement(IR_edge.sequence)
SSC_rc = reverse_complement(SSC_edge.sequence)

path1_seq = LSC_edge.sequence + IR_edge.sequence + SSC_edge.sequence + IR_rc
path2_seq = LSC_edge.sequence + IR_edge.sequence + SSC_rc             + IR_rc
```

**Biological rationale for path2:** The SSC region inverts independently
of the flanking IRs. `path2` captures the alternate SSC orientation:
the SSC is reverse-complemented while the flanking IR copies remain
identical to path1.

### 6.3 Sanity-check invariant

After construction, assert:

```
len(path1_seq) == len(path2_seq)
               == len(LSC) + 2 * len(IR) + len(SSC)
```

If this invariant fails, raise a `ValueError` — it indicates a
programming error in the construction step, not a bad input. Do not
catch this error; let it propagate so the process fails loudly.

### 6.4 FASTA output

```
write "{outdir}/path1.fasta":
    ">path1\n"
    path1_seq + "\n"

write "{outdir}/path2.fasta":
    ">path2\n"
    path2_seq + "\n"
```

No line-wrapping of the sequence. No other records in the file.

---

## 7. CLI arg-surface requirements

The script must be runnable from the command line for debugging, but
the **library entry point** (§8) is what C3 calls; the CLI is a
convenience wrapper.

Required CLI arguments:

| Argument | Form | Description |
|---|---|---|
| GFA file | Positional | Path to `assembly_graph.gfa` |
| `--outdir` | Optional flag, default `"."` | Directory for output FASTA files |
| `--json-out` | Optional flag | If given, write the Result metadata as JSON to this path |

**Do not** mirror the ptGAUL interface (`-e` / `-d` / `-o` triple).
Our implementation reads the GFA directly and does not accept
pre-extracted FASTA or pre-sorted depth files as inputs.

Example invocations:

```bash
# minimal
plastid_canonicalise.py assembly_graph.gfa

# with output
plastid_canonicalise.py assembly_graph.gfa --outdir results/ --json-out meta.json
```

The script should print a brief human-readable summary to stdout on
completion (e.g., `"branch=canonical lsc=edge_000001 ir=edge_000003 ssc=edge_000002"`).

---

## 8. Library entry point requirements

C3 calls C4 as a Python function, not as a subprocess. The entry
point must have the following signature:

```python
def canonicalise_plastid(
    gfa_path: str | Path,
    outdir: str | Path = ".",
) -> Result:
    ...
```

where `Result` is a namedtuple (or `dataclass`) with the fields
defined in §4.2.

**Contract:**

- `gfa_path` — path to `assembly_graph.gfa`; must exist and be
  readable. Raise `FileNotFoundError` if not.
- `outdir` — directory where `path1.fasta` and `path2.fasta` are
  written (only on `canonical` branch). The function must create this
  directory if absent (`os.makedirs(outdir, exist_ok=True)`).
- Return value — a `Result` namedtuple as specified in §4.2; never
  `None`.
- Side effects — writes FASTA files only on `canonical` branch; never
  on `resolved_circle` or `non_canonical`. Raises no exceptions for
  expected non-canonical inputs (edge-count mismatch, depth ties);
  raises `ValueError` only for genuinely corrupt inputs (duplicate
  edge names, `*` sequence, failed invariant).

The module should be importable without side effects:

```python
from plastid_canonicalise import canonicalise_plastid
result = canonicalise_plastid("assembly_graph.gfa", outdir="results/")
```

---

## 9. Test-case matrix

The tests in `scripts/tests/` (task 20 scope) must cover at minimum
the following cases. Each row shows the synthetic GFA to construct
and the assertions to make.

| # | Description | Synthetic GFA (S lines) | Assertions |
|---|---|---|---|
| 1 | 3-edge canonical | `edge_A len=90000 depth=1.0`; `edge_B len=25000 depth=2.0`; `edge_C len=15000 depth=1.0` | `result.branch == "canonical"`; `result.lsc_edge == "edge_A"`; `result.ir_edge == "edge_B"`; `result.ssc_edge == "edge_C"`; path1.fasta and path2.fasta exist |
| 2 | Uppercase `DP:f:` tag | Same as #1 but depth tag is `DP:f:2.0` on edge_B | `result.branch == "canonical"` (parser accepts uppercase) |
| 3 | Missing depth tag on IR | Same as #1 but edge_B has no depth tag → treated as depth 0.0; edge_A and edge_C also have no depth tags | `result.branch == "non_canonical"`; `result.non_canonical_reason` contains `"depth_tie"` |
| 4 | Missing depth tag on IR only | edge_A depth=1.0; edge_B has no depth tag (→ 0.0); edge_C depth=1.0 | `result.branch == "non_canonical"` because deepest edge is now A or C (tied at 1.0), and longer of those is A, so LSC=A, IR=A — LSC/IR collision triggers `"lsc_ir_collision"` |
| 5 | LSC–IR collision | edge_A len=90000 depth=3.0; edge_B len=25000 depth=2.0; edge_C len=15000 depth=1.0 | `result.branch == "non_canonical"`; reason is `"lsc_ir_collision"` (edge_A is both longest and deepest) |
| 6 | 1-edge resolved circle | Single `edge_A len=160000 depth=50.0` | `result.branch == "resolved_circle"`; no FASTA files written; `result.lsc_edge is None` |
| 7 | 2-edge | Two S lines | `result.branch == "non_canonical"`; `result.non_canonical_reason == "edge_count:2"` |
| 8 | 4-edge | Four S lines | `result.branch == "non_canonical"`; `result.non_canonical_reason == "edge_count:4"` |
| 9 | Round-trip length invariant | Same as #1 | `len(path1_seq) == len(path2_seq) == 90000 + 2*25000 + 15000 == 155000` |
| 10 | path2 uses rc(SSC) at correct offset | Same as #1 with known sequences | `path2_seq[90000+25000 : 90000+25000+15000] == reverse_complement(SSC_seq)` |

**Constructing synthetic GFAs for tests:**

Each test creates a minimal GFA string. Only S lines are required;
L lines can be omitted. Example for test #1:

```
H	VN:Z:1.0
S	edge_A	{LSC_seq_90k_chars}	dp:f:1.0
S	edge_B	{IR_seq_25k_chars}	dp:f:2.0
S	edge_C	{SSC_seq_15k_chars}	dp:f:1.0
```

Where `{LSC_seq_90k_chars}` is any DNA string of the stated length (e.g.,
`"A" * 90000`). Use distinct nucleotide content for each edge so that
sequence comparisons in tests #9 and #10 are meaningful.

---

## 10. Implementation-choice boundaries

The following details are intentionally left open; any choice within
the stated range is acceptable.

| Decision | Acceptable range |
|---|---|
| Sequence class | Plain Python `str`, Biopython `Seq`, or any object that supports `str()` and reverse-complement. If Biopython is used, `Bio.Seq.reverse_complement()` is preferred over a custom RC function. |
| Depth regex | Any case-insensitive regex that captures the float value from `dp:f:VALUE` or `DP:f:VALUE`, e.g. `re.search(r'[Dd][Pp]:f:([0-9.]+)', tag_string)`. The regex need not accept scientific notation; Flye does not emit it. |
| Tie-breaking | Any deterministic rule; the spec's suggestion (depth→length→lex) is a recommendation, not a requirement. The rule must be documented in a comment in the implementation. |
| Result type | `collections.namedtuple`, `typing.NamedTuple`, or `dataclass(frozen=True)`. |
| Empty-sequence handling | Raise `ValueError` or return `non_canonical` with reason `"zero_length_edge"`. Either is acceptable; the test matrix expects `non_canonical`. |
| Malformed depth float | If the depth tag value is present but not parseable as a float (e.g., `dp:f:nan`), treat as 0.0 or raise `ValueError`. Either; document the choice. |
| Line endings | Accept `\r\n` as well as `\n`; strip before splitting. |
| H / P / W / other GFA lines | Silently skip (do not parse, do not error). |

---

## 11. Provenance

The algorithm implemented by C4 was originally expressed in
`reference-material/ptgaul/combine_gfa.py`, part of the ptGAUL
pipeline (Xu et al. 2023; GitHub: `Bean061/ptGAUL`). The ptGAUL
repository ships no licence; default copyright ("all rights reserved")
applies to its source code. The algorithm itself — a bioinformatics
recipe identifying quadripartite regions by sequence length and read
depth — is not copyrightable expression.

`bin/plastid_canonicalise.py` is a **clean-room re-implementation**
written from this specification alone. The task 20 implementer does
not read `combine_gfa.py`. This separation is required by
[CONSTITUTION.md rule 12](../CONSTITUTION.md), which establishes a
preference for re-implementation over vendoring unlicensed upstream
code.

References:

- Original algorithm: `reference-material/ptgaul/combine_gfa.py`
  (retained for provenance; do not copy or import from this path)
- ptGAUL GitHub: `https://github.com/Bean061/ptGAUL`
- ptGAUL paper: Xu L, Dong Z, Fang L, et al. (2023). "ptGAUL: A
  pipeline for the assembly and classification of plant organellar
  genomes." *Molecular Ecology Resources*. (Cite as available; use
  GitHub URL as fallback if paper is not yet published at
  implementation time.)
- [CONSTITUTION.md rule 12](../CONSTITUTION.md)
