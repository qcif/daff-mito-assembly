# Task 72 — Migrate Nextflow code to the v26 strict syntax

**Status:** STUB — shelved on creation (2026-10-15). Flesh out before
executing.

## Overview

The workflow is written for, and pinned to, Nextflow 25.04.0 (CI and the
Azure launcher both set that version). From 25.10 onward Nextflow
parses config files and scripts with a stricter grammar (the "v2 syntax
parser"), and 26.04 makes it the default. Under 26.04 the workflow
cannot start: script compilation fails before any process runs.

A config-only fix was tried on 2026-10-15 and reverted (see "Config
changes" below). It only exposed the same problem in the scripts, and
committing half a migration would have left the code untested on either
version. This task does the whole move: config, scripts and launcher
together, so the pinned Nextflow
version can be moved forward deliberately rather than being stuck on
25.04. The goal is no change in behaviour: same channels, same
outputs and the same `metadata.json` / `run_manifest.json` content.
Only the source syntax changes.

## Config changes (tried 2026-10-15, then reverted)

With these changes, every profile passed `nextflow config` on 26.04.6.
Apply them again as part of this task.

1. **`nextflow.config`, `plugins {}`.** The strict parser rejects the
   `if (System.getenv('WF5_AZURE_POOL')) { id 'nf-azure…'; id
   'nf-cloudcache…' }` block ("Unexpected input: 'if'"). Plugins can't
   be declared per profile either. Fix: keep only `nf-schema` in the
   config, and have the launcher pass the full set on the command line:
   - add a constant such as
     `NXF_PLUGINS="nf-schema@…,nf-azure@…,nf-cloudcache@…"` beside
     `NXF_VERSION` in `deploy/azure/run-wf5-azure.sh`;
   - in `build_command()`, add `-plugins "${NXF_PLUGINS}"` straight
     after `run main.nf`.
   Update the comments in both files that say the pin lives in
   `nextflow.config`.
2. **`conf/azure.config`, env-var check.** The top-level
   `[...].each { if (!System.getenv(name)) throw ... }` loop is rejected
   ("Unexpected input: '['"). A top-level `if` and a closure inside an
   assignment were both tried, and both are rejected too. Fix: remove
   the loop and move the check to a `validateAzureEnv()` function in
   `main.nf`, called first in the entry workflow, before
   `validateParams()`. In pseudocode:

       if 'azure' not in workflow.profile.split(','): return
       missing = [v for v in REQUIRED if not env(v)]
       if missing: error "ERROR: -profile azure: environment variable(s)
                          not set: <missing> (run via the launcher; see
                          .env.sample and docs/azure/setup.md)"

   REQUIRED is the same seven names the config listed. Caveat: when
   tried on 25.04.0, Nextflow sets up storage and cloud-cache before any
   `main.nf` code runs. So a missing storage key, `WF5_AZURE_WORKDIR` or
   `WF5_AZURE_CACHE` fails before this check. nf-azure's own message for
   the storage key is clear, but the cloud-cache one isn't. The check
   only catches the Batch credentials and `WF5_AZURE_POOL`, and that
   path was never actually reached in testing, because fake credentials
   fail earlier. Before relying on it, decide whether the launcher's own
   checks (`REQUIRED_VARS` in `load_env()`) make it redundant.
3. **`main.nf`, string continuations.** Seven lines in the
   `validateParams()` error messages start with `+` on the continuation
   line inside `error ( ... )`. The strict parser rejects this
   ("Unexpected input: '+'"). Fix: put the `+` at the end of the
   previous line and keep the indentation as it is.

## Known breakages in scripts (from `nextflow lint main.nf modules subworkflows`)

With the config changes above applied: 24 errors in 2 files, plus 43
warnings in 3 files (not yet triaged):

- `main.nf`
  - `import groovy.json.JsonSlurper` is not allowed. Use the
    fully-qualified name inline, or the built-in JSON helpers.
  - `splitJson` is not defined.
  - The top-level `workflow.onComplete { ... }` handler can't sit among
    script declarations. Move it into the entry workflow's `onComplete:`
    section.
- `modules/local/collate.nf`
  - `optArg` is a closure defined in the `script:` block. The strict
    parser doesn't resolve it in the script string (19 uses). Hoist it
    to a script-level function, or build the arguments before the
    string.
- `subworkflows/local/annotation_scoring.nf`: the assignment to
  `annotation = JOIN_ANNOTATION_SCORES.out.annotation` is rejected
  (emit syntax).
- Warnings: `validation.monochromeLogs` /
  `validation.defaultIgnoreParams` are "unrecognised config option" on
  26.04. Check whether this is nf-schema 2.4.0 against the new config
  schema, or whether a newer nf-schema is needed.

## Open questions

1. **Target version.** Pin 26.04.x, or wait for 26.10 (already
   released)? Whichever is chosen must be set in all three places
   together: `.github/workflows/tests.yml`,
   `.github/workflows/integration.yml` and `NXF_VERSION` in
   `deploy/azure/run-wf5-azure.sh`. `manifest.nextflowVersion` in
   `nextflow.config` stays at a floor the code still works on.
2. **`nf-cloudcache`.** 26.04 reports "Plugin nf-cloudcache does not
   exist in any repository" for `nf-cloudcache@0.4.3`. Find out whether
   it is now bundled in core, renamed or re-versioned. The Azure
   `-resume` state depends on it (task 57).
3. **Support both versions during the move?** Should the code still run
   on 25.04 with `NXF_SYNTAX_PARSER=v2` while the version is being
   moved, so CI can switch in a separate commit from the code change?

## Constraints

- No change in behaviour or outputs (CONSTITUTION rule 18). Provenance
  fields in `metadata.json` / `run_manifest.json` must not change, apart
  from the Nextflow version they record (rule 16).
- Must keep working on the Azure Batch executor, not just locally
  (rule 17). The Azure launcher path needs a run before this task is
  done.
- No new test framework (rule 15). The existing three surfaces are
  enough to verify the move.

## Exit criteria (draft)

- `nextflow config -profile <p> .` parses for every profile
  (`stub`, `integration`, `docker`, `singularity`, and `azure` with its
  env vars set) on the target version.
- `nextflow lint main.nf modules subworkflows` reports zero errors on
  the target version, and its warnings have been triaged.
- On the target version, `-profile stub -stub-run` passes and
  `-profile integration` passes with unchanged biology assertions.
- CI runs on the target version and is green.
- At least one Azure launcher run (stub mode is enough) completes on the
  target version.
- `.nextflow.log` shows no deprecation warnings that come from our own
  code.
