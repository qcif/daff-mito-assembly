# Unit tests for wf5 custom logic

Each `bin/*.py` component has a corresponding test module here.
Tests run inside the `wf5-scripts` container image so the runtime
environment matches production exactly.

Mock external tool calls (`seqkit`, `seqtk`, `minimap2`, etc.) at the
subprocess boundary; test only the Python logic. Load the module
in-process via `importlib` (registering it in `sys.modules` so
`unittest.mock.patch()` can address it by name) rather than invoking
the script as a subprocess — a subprocess contributes zero measured
branch coverage under `scripts/pytest.sh`'s `coverage run`, which
traces only the parent process.

Run locally:
```bash
docker run --rm \
  -v $(git rev-parse --show-toplevel)/bin:/opt/wf5/bin:ro \
  -v $(git rev-parse --show-toplevel)/scripts/tests:/opt/wf5/tests:ro \
  -v $(git rev-parse --show-toplevel)/tests/fixtures:/opt/wf5/fixtures:ro \
  neoformit/daff-wf5-scripts:<tag> \
  python -m pytest /opt/wf5/tests -v
```

See plan.md §5a for the testing strategy.
