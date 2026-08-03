"""Unit tests for bin/parse_samplesheet.py (C1).

Each test loads parse_samplesheet.py in-process (importlib, module
registered in sys.modules so patch() can address it by name) rather
than invoking it as a subprocess — a subprocess contributes zero
measured branch coverage under `scripts/pytest.sh`, since `coverage
run` traces only the parent process (task 27). C1 is pure stdlib, so
no tool-boundary mocking is needed here, unlike coverage_gate.py.

`main()` raises SystemExit(1) (via `_fail()`) rather than returning an
error code. `run()` below catches that and folds it back into a
`returncode`/`stderr` pair (via pytest's `capsys` fixture) so the
existing subprocess-shaped assertions carry over unchanged.
"""

import importlib.util
import json
import sys
from pathlib import Path
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parent.parent.parent / 'bin' / \
    'parse_samplesheet.py'


def _load_module():
    spec = importlib.util.spec_from_file_location(
        'parse_samplesheet', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules['parse_samplesheet'] = module
    spec.loader.exec_module(module)
    return module


parse_samplesheet = _load_module()


class Result:
    """Mimics the subset of subprocess.CompletedProcess used below."""

    def __init__(self, returncode, stderr):
        self.returncode = returncode
        self.stderr = stderr


def run(samplesheet, data_dir, tmp_path, capsys):
    out = tmp_path / 'out.json'
    argv = [
        'parse_samplesheet.py',
        '--samplesheet', str(samplesheet),
        '--data-dir', str(data_dir),
        '--out', str(out),
    ]
    with patch.object(sys, 'argv', argv):
        try:
            parse_samplesheet.main()
            returncode = 0
        except SystemExit as exc:
            returncode = exc.code
    stderr = capsys.readouterr().err
    return Result(returncode, stderr), out


def make_data_dir(tmp_path, files):
    """Create empty files under tmp_path/data/ and return the data dir."""
    d = tmp_path / 'data'
    d.mkdir()
    for name in files:
        p = d / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.touch()
    return d


def write_csv(tmp_path, content):
    p = tmp_path / 'samples.csv'
    p.write_text(content)
    return p


# ---------------------------------------------------------------------------
# Happy-path tests
# ---------------------------------------------------------------------------

class TestValidSheet:
    def test_valid_multirow_sheet_emits_in_order(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['plant.fastq.gz', 'animal.fq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'PLANT-01,plant_pt,plant.fastq.gz\n'
            'ANIMAL-01,animal_mt,animal.fq.gz\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr
        rows = json.loads(out.read_text())
        assert len(rows) == 2
        assert rows[0]['sample_id'] == 'PLANT-01'
        assert rows[0]['assembly_target'] == 'plant_pt'
        assert rows[1]['sample_id'] == 'ANIMAL-01'
        assert rows[1]['assembly_target'] == 'animal_mt'

    def test_valid_multifile_row_resolves_all_reads(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz', 'b.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'MULTI-01,animal_mt,a.fastq.gz|b.fastq.gz\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr
        rows = json.loads(out.read_text())
        assert len(rows[0]['reads']) == 2

    def test_reads_resolved_to_absolute_paths(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,plant_pt,a.fastq.gz\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr
        rows = json.loads(out.read_text())
        assert Path(rows[0]['reads'][0]).is_absolute()

    def test_optional_columns_absent_from_header_ok(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,plant_pt,a.fastq.gz\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr
        row = json.loads(out.read_text())[0]
        assert row['sample_info'] == ''
        assert row['sample_type'] == ''
        assert row['storage_location'] == ''

    def test_optional_columns_empty_cells_ok(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads,sample_info,sample_type,'
            'sample_receipt_date,storage_location\n'
            'S1,plant_pt,a.fastq.gz,,,,\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr

    def test_assembly_target_case_normalised_to_lowercase(
            self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,Plant_PT,a.fastq.gz\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr
        assert json.loads(out.read_text())[0]['assembly_target'] == 'plant_pt'

    def test_bad_date_warns_but_passes_through(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads,sample_receipt_date\n'
            'S1,plant_pt,a.fastq.gz,24/07/2026\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr
        assert 'WARNING' in result.stderr
        assert json.loads(out.read_text())[0]['sample_receipt_date'] == \
            '24/07/2026'

    def test_valid_iso_date_no_warning(self, tmp_path, capsys):
        """The success arm of date.fromisoformat — no warning emitted."""
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads,sample_receipt_date\n'
            'S1,plant_pt,a.fastq.gz,2026-07-24\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr
        assert 'WARNING' not in result.stderr
        assert json.loads(out.read_text())[0]['sample_receipt_date'] == \
            '2026-07-24'

    def test_all_extensions_accepted(self, tmp_path, capsys):
        exts = ['a.fastq', 'b.fastq.gz', 'c.fq', 'd.fq.gz']
        d = make_data_dir(tmp_path, exts)
        reads = '|'.join(exts)
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            f'S1,plant_pt,{reads}\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr
        assert len(json.loads(out.read_text())[0]['reads']) == 4

    def test_all_three_assembly_targets_accepted(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,animal_mt,a.fastq.gz\n'
            'S2,plant_pt,a.fastq.gz\n'
            'S3,plant_mt,a.fastq.gz\n'
        ))
        result, out = run(csv, d, tmp_path, capsys)
        assert result.returncode == 0, result.stderr
        targets = [r['assembly_target'] for r in json.loads(out.read_text())]
        assert targets == ['animal_mt', 'plant_pt', 'plant_mt']


# ---------------------------------------------------------------------------
# Failure tests
# ---------------------------------------------------------------------------

class TestFailures:
    def test_duplicate_sample_id_fails(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'DUP-01,plant_pt,a.fastq.gz\n'
            'DUP-01,animal_mt,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'DUP-01' in result.stderr

    def test_absolute_reads_path_fails(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, [])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,plant_pt,/absolute/path.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'absolute path' in result.stderr

    def test_missing_reads_file_fails(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, [])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,plant_pt,nonexistent.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'not found' in result.stderr

    def test_bad_assembly_target_fails(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,fungal_mt,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'fungal_mt' in result.stderr

    def test_bare_kingdom_value_rejected(self, tmp_path, capsys):
        """A bare kingdom name (no organelle) is not a valid target."""
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,plant,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'plant' in result.stderr

    def test_unknown_column_in_header_fails(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads,extra_column\n'
            'S1,plant_pt,a.fastq.gz,value\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'extra_column' in result.stderr

    def test_missing_required_column_fails(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, [])
        csv = write_csv(tmp_path, (
            'sample_id,reads\n'
            'S1,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'assembly_target' in result.stderr

    def test_legacy_kingdom_column_rejected(self, tmp_path, capsys):
        """A samplesheet using the pre-migration 'kingdom' column must be
        rejected — both because assembly_target is missing (required) and
        because 'kingdom' is now an unknown column."""
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'S1,plant,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        # Both errors surface (all errors collected before exit)
        assert 'assembly_target' in result.stderr
        assert 'kingdom' in result.stderr

    def test_bad_sample_id_pattern_fails(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'bad id!,plant_pt,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0

    def test_unsupported_reads_extension_fails(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.bam'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,plant_pt,a.bam\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'extension' in result.stderr

    def test_all_errors_collected_before_exit(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        # 3 distinct problems: bad target, duplicate id, absolute path
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,fungal_mt,a.fastq.gz\n'
            'S1,plant_pt,/abs/path.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert result.stderr.count('ERROR:') >= 3

    def test_no_header_row_fails(self, tmp_path, capsys):
        """An empty samplesheet has no fieldnames at all."""
        d = make_data_dir(tmp_path, [])
        csv = write_csv(tmp_path, '')
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'empty' in result.stderr.lower()

    def test_empty_reads_cell_fails(self, tmp_path, capsys):
        d = make_data_dir(tmp_path, [])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,plant_pt,\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'reads column is empty' in result.stderr

    def test_stray_pipe_in_reads_fails(self, tmp_path, capsys):
        """A double pipe leaves an empty entry after split('|')."""
        d = make_data_dir(tmp_path, ['a.fastq.gz', 'b.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,assembly_target,reads\n'
            'S1,plant_pt,a.fastq.gz||b.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path, capsys)
        assert result.returncode != 0
        assert 'stray pipes' in result.stderr
