---
name: azure-run
description: Run the workflow on Azure Batch against real samples given as az:// blob URIs
---

Run the pipeline on Azure Batch for "$1": one or more samples whose reads
are already in blob storage, each given as an `az://<container>/<path>`
URI. Nothing is uploaded from this machine. Nextflow runs locally as the
launcher and each process runs as a Batch task on a named pool.

Background: docs/azure/README.md (layout, pools, debugging) and
tasks/completed/57_azure_batch_deployment.md (what was verified, and the
fixes it needed). Read those before improvising anything not covered here.

## Facts (verified in task 57; re-check rather than assume if something fails)

- **Always launch through `deploy/azure/run-wf5-azure.sh`.** It sources
  `.env.azure`, checks the pool and storage first, sets the per-run
  `WF5_AZURE_*` vars `conf/azure.config` needs, and pins `NXF_VER=25.04.0`.
  Never call `nextflow run -profile azure` by hand.
- **Credentials live in `.env.azure` (git-ignored).** Never print it, never
  echo a key, and never put a key on a command line that you show in
  output. Source it inside the command, as in the snippets below.
- **Pools:** `wf5_test` is the default. Use `wf5` (production) only when the
  user asks for it (`--pool wf5`). Both are `Standard_L8as_v3` (8 vCPU / 64
  GiB), so they need no resource clamp like `remote_resources.config`.
- **Refs** are staged per task from `az://refdata-wf5/<version>/`, and
  `--refs-version` is required. A version must exist as an extracted
  bundle on blob, not only as a local `refs/<version>/`.
- **Work dir:** `az://workdata/wf5/<run id>/`, one per launch, removed by
  the lifecycle rule 14 days after the last write. Cache:
  `az://cache/wf5/{test,prod}/`.
- **Outputs are published locally**, to the launcher's `--outdir`. Nextflow
  copies them down from blob.
- **Expected duration:** the 3-sample integration run took about 43 min.
  `METAFLYE` alone can take up to about 20 min per sample, and container
  pulls and Batch scheduling add more. Node scale-up adds a few minutes
  before the first task starts.

## Procedure

1. **Check the URIs exist.** For each `az://<container>/<blob>`, check it
   before writing a samplesheet. `main.nf` checks again at intake, but a
   typo found here costs seconds, not a pool spin-up:
   ```
   ( set -a
     . ./.env.azure
     set +a
     az storage blob show \
       --account-name "$STORAGE_ACCOUNT_STD" \
       --account-key "$AZURE_STORAGE_ACCOUNT_KEY" \
       --container-name <container> --name <blob> \
       --query "{name: name, size: properties.contentLength}" -o tsv )
   ```
   If the user gives a prefix rather than files, list it with
   `az storage blob list ... --prefix <path>/` and confirm with the user
   which files belong to which sample. Don't guess the grouping.

2. **Samplesheet.** Build one per spec/00-overview.md §0. The required
   columns are `sample_id`, `assembly_target` and `reads`. The optional
   ones are `sample_info`, `sample_type`, `sample_receipt_date` and
   `storage_location`. An unknown column fails the parse.
   - `reads` holds the full `az://` URI. Separate multiple files for one
     sample with `|`. A sheet is either all-URI or all-relative, never
     mixed, and only the `az://` scheme is accepted.
   - Save it as `tests/manual/samples_<run-name>.csv`, beside the
     existing sheets.
   - Confirm `assembly_target` for each row with the user, per the
     CONSTITUTION's hard constraint 2. It is declared at intake and never
     inferred.
   - **Sample ID = the blob's filename** (the last path segment), minus its
     FASTQ extension (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`):
     `.../barcode21_SAC1451-2-1.fastq.gz` → `barcode21_SAC1451-2-1`. Never
     shorten it to `barcodeNN`. Barcodes are reused across runs, so
     `barcode21` alone is ambiguous. Sanitise only as far as `sample_id`
     requires (`[A-Za-z0-9_.-]+`, spec/00 §0): drop browser suffixes like
     ` (1)` and replace any other disallowed character with `_`. For a
     multi-file sample, use the name the files share, and ask if they
     share none.

3. **Pick the refs version.** By default, use the newest version published
   to blob:
   ```
   ( set -a
     . ./.env.azure
     set +a
     az storage blob list \
       --account-name "$STORAGE_ACCOUNT_STD" \
       --account-key "$AZURE_STORAGE_ACCOUNT_KEY" \
       --container-name refdata-wf5 --delimiter / \
       --query "[].name" -o tsv )
   ```
   A version with no `<version>/manifest.json` was never fully published,
   because `scripts/publish_refs_extracted.sh` uploads the manifest last.
   Don't use it.

4. **Stub check, when something structural changed.** It catches a wiring,
   auth or pool problem before a long run. It takes about 16 min on Azure,
   so skip it when only `bin/` or report code changed since the last
   good Azure run:
   ```
   deploy/azure/run-wf5-azure.sh --stub --refs-version <ver> --yes \
     --outdir tests/manual/output/stub_check_<run-name>
   ```
   One `TRANSLATE_ANNOTATION_CDS` ignore here is a known issue (task 62,
   an empty directory output), not a new failure. Delete the stub outdir
   afterwards.

5. **Launch.** The launcher runs Nextflow in the foreground on this
   machine, so start it with `run_in_background` and `--yes`, and send its
   output to a per-run log:
   ```
   deploy/azure/run-wf5-azure.sh \
     --refs-version <ver> \
     --samplesheet tests/manual/samples_<run-name>.csv \
     --outdir tests/manual/output/<run-name> \
     --yes \
     > tests/manual/output/run_<run-name>.log 2>&1
   ```
   - **No `--data-dir`.** With an all-URI sheet, the launcher reports "Data
     dir: (none; sheet lists az:// URIs)". If it reports a byte count to
     upload, stop: the sheet isn't all-URI.
   - Pass `--pool wf5` or `--polish` only when asked. Send extra pipeline
     params after `--`. For example, for the report's run-context panel:
     `-- --facility "<lab>" --analyst_name "<name>"`. `.env.azure` may
     define `FACILITY_NAME` / `ANALYST_NAME`, but the launcher does not
     pass them on.
   - **One Azure run at a time from this checkout.** The launcher runs from
     the repo root and writes `.nextflow.log` and `.nextflow/` there, so
     two concurrent runs would overwrite each other's log and history.
   - The per-run outdir must be new. Check `ls tests/manual/output/` first.

6. **Wait on it properly.** Use a Monitor with an until-loop on the log.
   Don't poll or sleep:
   ```
   until grep -qE 'Workflow completed|Workflow failed|ERROR( ~|:)' \
     tests/manual/output/run_<run-name>.log 2>/dev/null \
     || ! pgrep -f '[r]un-wf5-azure.sh' >/dev/null; do sleep 30; done
   tail -c 3000 tests/manual/output/run_<run-name>.log | tr '\r' '\n' | tail -40
   ```
   Match the failure markers and the launcher's own exit (the `pgrep` arm),
   not just success. Silence must never be mistaken for progress. The
   `[r]` stops `pgrep -f` matching the Monitor's own command line, which
   contains the same string. A
   Monitor expires after at most 30 min, so expect to re-arm at least once.
   When it expires, tail the log to confirm tasks are still moving before
   re-arming. If the log shows no submissions for a long time, check the
   pool with `source deploy/azure/batch-helpers.sh && az_pool_show` (node
   counts, allocation errors).

7. **Keep the session log.** Copy `.nextflow.log` to
   `tests/manual/output/<run-name>/nextflow.log` straight away, before
   any later run rotates it. It holds the `Session UUID` that `-resume`
   needs (see Resuming) and the work-dir hashes that `az_fetch_nf_logs`
   reads.

8. **Review.** The outputs are already in `tests/manual/output/<run-name>/`
   (git-ignored). Report each sample's terminal status (`ok`,
   `low_coverage`, `no_assembly`, `no_barcode`, `fail` or `error`) and
   keep them distinct, per CONSTITUTION principle 7. `no_assembly`
   (assembled nothing, or binned off-target) and `no_barcode` (assembled
   fine, but no locus could be extracted) are different findings. `fail`
   (a deliberate coverage-gate decision) is not `error` (the pipeline
   broke).

   Status is in each `<sample>/metadata.json` as `sample_status` and
   `sample_status_reason`. The gate numbers are under `coverage.gate`, the
   assembly summary under `assembly`, and the loci under `barcodes.loci`.
   Skip `assembly`'s `contigs` and `bin_metadata`, which are huge. Print
   the fields you need. Dumping whole sections floods the output.

   **Check that every sheet row has a `metadata.json`.** A missing one
   means the sample crashed or was dropped. On Azure, the known cause is an
   empty directory output (task 62). Report it as such and don't fold it
   into `no_assembly`.

   Point the user at `run-report.html` and `run_manifest.json` first. The
   per-sample `report.html` comes from `COLLATE`.

9. **Clean up.** Trim stale `.nextflow.log.N` in the repo root (keep
   `.nextflow.log` and `.1` / `.2`). Leave the blob work dir alone: it
   expires on its own, and `-resume` needs it until then. Delete
   `workdata/wf5/<run id>/` only if the user asks. Never touch other
   workflows' prefixes in the shared storage account (spec/04a).

## If a task fails

Per-sample crash isolation is on (`errorStrategy = 'ignore'`), so one
ignored process is often a real biological outcome rather than a bug.
To see what happened, pull the task logs:
```
( set -a
  . ./.env.azure
  set +a
  source deploy/azure/batch-helpers.sh
  az_fetch_nf_logs tests/manual/output/<run-name>/nextflow.log \
    tests/manual/output/<run-name>/task_logs )
```
That prints the tasks with a non-empty `.command.err` and the non-zero
exit codes. `az_job_logs` shows Batch-side task states. A failure that
only happens on Azure, and not locally, is often a staging problem: an
input not passed through a channel (rule 13, spec/01 §1a), or an empty
directory output (task 62). Don't patch pipeline code mid-run. Document
the failure and raise a task.

## Resuming

`-resume` is cheap after a report-side or config fix, because assembly
stays cached. It reaches back only as far as the work dir's 14-day
lifecycle, and ignored-error tasks re-run anyway. It needs the original
session id (the `Session UUID` in the saved `nextflow.log`). If
`.nextflow/history` no longer has that session, seed it as in
docs/azure/README.md "Resuming a run", then:
```
deploy/azure/run-wf5-azure.sh <same args> -- -resume <session id>
```
Don't rely on it to pick up in-place edits to `bin/*.py` (see
tasks/todo.md). Run clean after Python changes.
