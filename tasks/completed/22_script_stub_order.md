# Task 22 — Reorder `script:` / `stub:` blocks in every NF module

**Phase:** Cross-cutting hygiene (no phase; touches every module).
**Goal:** Clear the
`Invalid process definition -- check for missing or out-of-order section labels`
error that the nf-language-server (VS Code extension) raises on every
process in [`modules/local/`](../modules/local/).

Root cause: all 18 modules currently place `stub:` **before** `script:`.
The Nextflow language specification requires the main execution block
(`script:` / `shell:` / `exec:`) to precede the `stub:` block. The
Nextflow runtime is permissive enough that `-stub-run` still works,
which is how the codebase drifted into the wrong order — but the
language server (and any downstream `nextflow lint` tooling) is strict.

Canonical intra-process order per the Nextflow spec:

1. Directives (`tag`, `label`, `publishDir`, `errorStrategy`, …)
2. `input:`
3. `output:`
4. `when:` (optional)
5. `script:` / `shell:` / `exec:`   ← main block, comes first
6. `stub:`                          ← comes AFTER `script:`

**Exit criteria:**

- Every module in
  [`modules/local/`](../modules/local/) has `script:` on a lower line
  number than `stub:`.
- Opening any module in VS Code with the Nextflow extension raises
  no "out-of-order section labels" error.
- `nextflow run . -profile stub,docker -stub-run` still hits 17/17
  process stubs successfully — no change in observable pipeline
  behaviour.
- No changes to `main.nf`, `nextflow.config`, `conf/*.config`,
  process directives, `input:`/`output:` blocks, or the command
  bodies inside `script:`/`stub:`.

**Not in scope:**

- No refactor of directives, resource labels, or publish semantics.
- No changes to the command bodies (bash inside `script:` / `stub:`
  is copied byte-for-byte, only its position moves).
- No conversion of any module to `shell:` or `exec:` — modules that
  currently use `script:` stay on `script:`.
- No migration to nf-core module conventions.

**Cross-cutting rules (from
[spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints) and
[CONSTITUTION.md](../CONSTITUTION.md)):**

- The fix is behaviourally inert; per
  [CONSTITUTION rule 13](../CONSTITUTION.md) (small, reversible changes
  land as isolated commits) this ships as a single mechanical commit
  separate from any feature work.
- Groovy preludes (`def foo = params.bar[meta.assembly_target]`) that
  currently sit at the top of `script:` stay inside the moved
  `script:` block — the block relocates as a unit, `def` line first.

---

## 1. Affected files (18 total)

All modules under [`modules/local/`](../modules/local/):

`annotate.nf`, `bandage_ng.nf`, `bin_target.nf`, `blast_validate.nf`,
`chopper.nf`, `collate.nf`, `coverage_gate.nf`, `filtlong.nf`,
`medaka.nf`, `metaflye.nf`, `miniprot_extract.nf`,
`nanoplot_clean.nf`, `nanoplot_raw.nf`, `organelle_map.nf`,
`parse_samplesheet.nf`, `recruit.nf`, `run_report.nf`,
`validate_samplesheet.nf`.

Audit command to confirm the current state before starting:

```bash
grep -n -E '^\s*(script|stub):' modules/local/*.nf
```

After the fix, every file's `script:` line number must be lower than
its `stub:` line number.

## 2. Mechanical fix

For each module, swap the two blocks so `script:` appears first and
`stub:` follows. Preserve all whitespace, comments, and the exact
text inside each block.

### Example — before ([`bandage_ng.nf`](../modules/local/bandage_ng.nf))

```groovy
    output:
    tuple val(meta), path(assembly), path(gfa), path(info),
          path("${meta.sample_id}.graph.png"), emit: assembly

    stub:
    """
    touch ${meta.sample_id}.graph.png
    """

    script:
    """
    BandageNG image ${gfa} ${meta.sample_id}.graph.png \\
        --height 600 \\
        --width 800
    """
}
```

### Example — after

```groovy
    output:
    tuple val(meta), path(assembly), path(gfa), path(info),
          path("${meta.sample_id}.graph.png"), emit: assembly

    script:
    """
    BandageNG image ${gfa} ${meta.sample_id}.graph.png \\
        --height 600 \\
        --width 800
    """

    stub:
    """
    touch ${meta.sample_id}.graph.png
    """
}
```

Note: modules with a Groovy prelude keep the `def` inside `script:`
when it moves. Two examples in the current tree:

- [`coverage_gate.nf`](../modules/local/coverage_gate.nf) —
  `def limits = params.coverage_limits[meta.assembly_target]` stays
  as the first line of the moved `script:` block.
- [`bin_target.nf`](../modules/local/bin_target.nf) —
  `def codes = params.genetic_code_tables[meta.assembly_target].join(',')`
  stays as the first line of the moved `script:` block.

## 3. Verification

1. **IDE lint** — reopen each module in VS Code. The
   "Invalid process definition" error must be gone on all 18 files.
2. **Grep audit** — the following command must return no lines
   (i.e. every module has `script:` before `stub:`):
   ```bash
   awk '
     /^\s*script:/  { script=NR }
     /^\s*stub:/    { stub=NR }
     ENDFILE {
       if (script && stub && stub < script) print FILENAME
       script=0; stub=0
     }' modules/local/*.nf
   ```
3. **Stub run** — behaviour must be unchanged:
   ```bash
   nextflow run . -profile stub,docker -stub-run
   ```
   Expect 17/17 processes to hit their stub and complete, matching
   the pre-fix output.
4. **Integration run (optional, belt-and-braces)** — only if
   convenient, since the change is behaviourally inert:
   ```bash
   nextflow run . -profile integration
   bash tests/integration/assertions.sh
   ```
   Expect no change vs. the current baseline.

## 4. Test-fixture / unit-test impact

**None.** No integration fixtures change (no process output changes).
No Python code is touched, so no unit test changes are needed. The
existing integration
[`assertions.sh`](../tests/integration/assertions.sh) suite is the
regression check for this task — it must continue to pass unchanged.

## 5. Deliverables checklist

- [x] All 18 modules under
      [`modules/local/`](../modules/local/) have `script:` before
      `stub:` in file order.
- [x] Grep audit (§3 command) returns no results.
- [ ] VS Code Nextflow extension reports no "out-of-order section
      labels" errors on any module.
- [ ] `nextflow run . -profile stub,docker -stub-run` — 17/17
      processes complete, no regressions.
- [x] Single mechanical commit with message along the lines of
      *"Reorder script/stub blocks for nf-language-server compliance"*,
      not bundled with feature work.

## 6. Outcomes

- Files touched: 18 / 18
- Stub-run result: pending
- IDE lint result: pending (requires VS Code reload)
- Any incidental deviation (e.g. a module with a genuinely different
  structural issue surfaced by the reordering): none — all modules
  had a straightforward stub/script pair to swap.
