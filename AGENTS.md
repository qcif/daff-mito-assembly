# Organelle genome assembly workflow

- Read the CONSTITUTION.md and brief.md
- plan.md is the root document for the project specification, which lives in
  ./spec/ and should be consulted for all design decisions. At the very least,
  read spec/00-overview.md and spec/01-pipeline-flow.md.
- Task definitions live in ./tasks/, and completed tasks can be seen in
  ./tasks/completed/.
- Workflow development should only be done with the guidance of an accompanying
  ./task/*.md brief.
- Run `flake8` on `bin/*.py` / `scripts/tests/*.py` using the `claude`
  venv, not Docker — activate with
  `. /home/cameron/.local/envs/claude/bin/activate`. flake8 only
  parses the source, so it doesn't need the pinned runtime deps.
- To run `pytest` / `coverage` for `bin/*.py` unit tests: `scripts/pytest.sh`
