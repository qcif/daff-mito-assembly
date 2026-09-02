"""Result objects for the report template.

Two reusable patterns:

- `AbstractDataRow`: one dict row cast to typed attributes, driven by a
  static `COLUMNS` list on the subclass. Use for fixed, single-row things
  like sample metadata or a per-sample QC summary.

- `AbstractResultRows`: many dict rows cast per column type, driven by a
  `COLUMN_METADATA` dict loaded from a CSV schema in `schema/`. Use for
  DataTables-style tabular results (e.g. barcode extraction hits).

Subclasses below are stubs — flesh out COLUMNS / COLUMN_METADATA and add
any result-specific helpers (bootstrap classes, chart-data derivations,
etc.) as pipeline stages are implemented.
"""

import csv
from typing import Union, get_args, get_origin


class FLAGS:
    SUCCESS = 'success'
    WARNING = 'warning'
    DANGER = 'danger'
    NONE = 'secondary'


def colour_to_bs_class(colour: str) -> str:
    """Convert a colour string to a Bootstrap contextual class."""
    mapping = {
        'green': FLAGS.SUCCESS,
        'yellow': FLAGS.WARNING,
        'red': FLAGS.DANGER,
    }
    return mapping.get(colour.lower(), FLAGS.NONE)


def _csv_to_dict(csv_path, index_col='colname'):
    """Load a schema CSV into an ordered dict keyed by `index_col`."""
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)
        ordered_colnames = [row[0].strip() for row in reader]
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        data = {
            row[index_col].strip(): dict(row.items())
            for row in reader
        }
    return {name: data[name] for name in ordered_colnames}


class AbstractDataRow:
    """Single-row result cast to typed attributes.

    Subclasses declare `COLUMNS = [(name, type), ...]`. Values missing from
    the source row are set to None; `Optional[T]` unwraps to `T`.
    """

    COLUMNS: list[tuple[str, type]] = []

    def __init__(self, row: dict):
        for colname, _type in self.COLUMNS:
            raw = row.get(colname)
            if raw is None:
                value = None
            else:
                origin = get_origin(_type)
                if origin is Union:
                    allowed = [
                        t for t in get_args(_type) if t is not type(None)
                    ]
                    value = allowed[0](raw.strip()) if allowed else None
                else:
                    value = _type(raw.strip())
            setattr(self, colname, value)

    def to_json(self):
        return {name: getattr(self, name) for name, _ in self.COLUMNS}


class AbstractResultRows:
    """Tabular result driven by a `COLUMN_METADATA` schema dict.

    COLUMN_METADATA rows may carry:
      - `type`: one of int / float / scientific / bool (used for casting)
      - `label`: display label; empty means the column is hidden
      - `primary_display`: truthy means show in the compact/default view
      - anything else the subclass wants to consult

    Subclasses may override `__init__` to derive extra per-row fields
    (bootstrap class, null-row handling, etc).
    """

    COLUMNS: list[str] = []
    COLUMN_METADATA: dict = {}

    def __init__(self, rows):
        self.rows = self._parse_rows(rows)
        self.columns_display = [
            c for c in self.COLUMN_METADATA
            if self.COLUMN_METADATA[c].get('label')
        ]
        self.columns_primary_display = [
            c for c in self.columns_display
            if self.COLUMN_METADATA[c].get('primary_display')
        ]

    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def __getitem__(self, index):
        return self.rows[index]

    @classmethod
    def from_csv(cls, path, delimiter='\t'):
        with path.open() as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            return cls(list(reader))

    def _parse_rows(self, rows):
        return [
            {
                colname: self._cast(
                    row[colname].strip() if row[colname] is not None else '',
                    self.COLUMN_METADATA.get(colname, {}).get('type'),
                )
                for colname in self.COLUMNS
                if colname in row
            }
            for row in rows
        ]

    def _cast(self, value, type_str):
        try:
            if not type_str:
                return value
            if type_str == 'int':
                return f"{int(float(value)):,}"
            if type_str == 'float':
                if value.lower() in ('na', 'nan', 'n/a'):
                    return None
                return float(value)
            if type_str == 'scientific' and 'e' in value:
                return f"{float(value):.2e}"
            if type_str == 'bool':
                return value.lower() in ('true', '1', 'yes', 'y')
        except ValueError:
            return None
        return value

    def to_json(self):
        return self.rows


# ---------------------------------------------------------------------------
# `metadata.json` (task 42) is the renderer's only input — `sample_id`,
# `kingdom`, `sample_status` and the submitter-supplied optional columns
# come straight from the parsed dict (report.py's `build_context`), so
# there is no per-stage `AbstractDataRow` subclass here. `AbstractDataRow`
# / `AbstractResultRows` stay as the reusable base classes 43b's tabular
# stage results (barcode hits, per-contig assembly stats, ...) build on.
#
# Example tabular-result stub. Delete or clone as needed.
#
# class BarcodeHits(AbstractResultRows):
#     COLUMN_METADATA = _csv_to_dict(SCHEMA.EXAMPLE_FIELD_CSV)
#     COLUMNS = list(COLUMN_METADATA.keys())
