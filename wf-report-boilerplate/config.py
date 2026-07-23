"""Configuration for the report module.

Result-file locations are resolved by glob against a per-sample result
directory. Add new `@property` accessors as new stages produce outputs.
"""

import csv
import os
from datetime import datetime
from functools import cached_property
from pathlib import Path

import yaml

from .utils import path_safe

PARENT_DIR = Path(__file__).parent
ROOT_DIR = Path(__file__).parents[2].resolve()

# TODO: update once the wf5 repo URL exists.
REPO_URL = 'https://github.com/qcif/daff-biosecurity-wf5'
WF_NAME = 'WF5 (Organelle Assembly + Barcode Recovery)'


class Config:

    TIMESTAMP_FILE = '*_start_timestamp.txt'
    VERSIONS_PATH = ROOT_DIR / 'versions.yml'
    DEFAULT_PARAMS_PATH = ROOT_DIR / 'params/default_params.yml'
    FLAGS_CSV = ROOT_DIR / 'flags.csv'

    class SCHEMA:
        # Column-metadata CSVs describing tabular result schemas. Add one
        # per tabular result type; consumed by AbstractResultRows subclasses.
        EXAMPLE_FIELD_CSV = PARENT_DIR / 'schema/example_fields.csv'

    class OUTPUTS:
        REPORT_FILE_TEMPLATE = '{sample_id}_report.html'

    class REPORT:
        TITLE = "Organelle assembly + barcode recovery report"
        SUBTITLE = (
            "Results generated through the"
            f' <a href="{REPO_URL}" target="_blank">'
            f'{WF_NAME} Nextflow pipeline</a>.')

    @property
    def result_dir(self) -> Path:
        """Ensure that result dir is propagated between Config instances."""
        return Path(os.environ['RESULT_DIR'])

    @property
    def report_path(self) -> Path:
        return (self.result_dir / self.OUTPUTS.REPORT_FILE_TEMPLATE.format(
            sample_id=path_safe(self.sample_id),
        )).absolute()

    # ---- Result file accessors -------------------------------------------
    # Add one @property per output artefact the report needs to read.
    # Use `_get_file_by_pattern` with the glob emitted by the workflow.
    #
    # Example:
    #
    # @property
    # def assembly_fasta_path(self) -> Path:
    #     return self._get_file_by_pattern('*_organelle.fasta')
    # ----------------------------------------------------------------------

    @cached_property
    def sample_id(self) -> str:
        """Derive the sample ID from a canonical result file.

        Stub: replace the pattern with a wf5 result file (e.g. the assembly
        FASTA or a per-sample manifest) once one is defined.
        """
        # TODO: pick a canonical wf5 file and derive sample_id from its name.
        raise NotImplementedError(
            "sample_id derivation not yet wired for wf5 outputs"
        )

    @cached_property
    def start_time(self):
        """Return workflow start timestamp, or None if unavailable."""
        path = self._get_file_by_pattern(
            self.TIMESTAMP_FILE, ignore_missing=True,
        )
        if path and path.exists():
            return datetime.strptime(path.read_text().strip(), '%Y%m%d%H%M%S')
        return None

    @cached_property
    def flags(self) -> list[dict[str, str]]:
        """Load flag definitions from a CSV file, if present."""
        if not self.FLAGS_CSV.exists():
            return []
        with self.FLAGS_CSV.open() as f:
            return list(csv.DictReader(f))

    @cached_property
    def default_params(self) -> dict:
        if not self.DEFAULT_PARAMS_PATH.exists():
            return {}
        with self.DEFAULT_PARAMS_PATH.open() as f:
            return yaml.safe_load(f) or {}

    def load(self, result_dir: Path):
        """Bind the result directory used by all subsequent lookups."""
        os.environ['RESULT_DIR'] = str(result_dir)

    def _get_file_by_pattern(
        self,
        file_pattern: str,
        ignore_missing: bool = False,
    ) -> Path:
        paths = list(self.result_dir.glob(file_pattern, case_sensitive=False))
        if paths:
            return paths[0]
        if ignore_missing:
            return None
        raise FileNotFoundError(
            f'No file matching pattern: {self.result_dir / file_pattern}'
        )
