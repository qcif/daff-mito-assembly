# tests/unit/

Unit tests for the custom-logic components in `bin/`. One test module per
script, using Python `unittest`. Tests run inside the `wf5-scripts` container
so the runtime environment matches production.

See plan.md §5a for the full testing strategy.

## Run locally

```bash
python -m unittest discover -s tests/unit -p 'test_*.py'
```

## Adding tests

Create `test_<component>.py` alongside each new `bin/<component>.py`.
Mock external tool calls (seqkit, seqtk, minimap2, etc.) at the boundary;
test only the Python logic.
