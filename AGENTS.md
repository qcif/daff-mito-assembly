# Organelle genome assembly workflow

- Read the CONSTITUTION.md and brief.md
- plan.md is the root document for the project specification, which lives in
  ./spec/ and should be consulted for all design decisions. At the very least,
  read spec/00-overview.md and spec/01-pipeline-flow.md.
- Task definitions live in ./tasks/, and completed tasks can be seen in
  ./tasks/completed/.
- Workflow development should only be done with the guidance of an accompanying
  ./task/*.md brief.

## Linting and testing
- Run `flake8` on `bin/*.py` / `scripts/tests/*.py` using the `claude`
  venv, not Docker — use it directly with
  `/home/cameron/.local/envs/claude/bin/flake8`, no need to activate the venv.
- To run `pytest` / `coverage` for `bin/*.py` unit tests: `scripts/pytest.sh`.
  This script should be your first choice, but you can also use Claude's python
  binary directly if required:
  `/home/cameron/.local/envs/claude/bin/python -m pytest ...`

## Running Nextflow

- Remove any .nextflow.log.N files where N>2
- If using -resume, make sure ./work/ only contains data you want to re-use.
  Delete the whole dir when you want a clean run to stop stale task data from
  accumulating.
- To clean up after a run: `.claude/scripts/clean_nextflow_run.sh`. Don't use a
  Bash `rm` for this, use the standard process.
