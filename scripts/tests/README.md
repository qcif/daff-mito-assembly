# Unit tests for wf5 custom logic

Each `bin/*.py` component has a corresponding test module here.
Tests run inside the `wf5-scripts` container image so the runtime
environment matches production exactly.

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
