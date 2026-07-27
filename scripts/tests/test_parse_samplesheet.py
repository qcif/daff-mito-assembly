"""Unit tests for bin/parse_samplesheet.py (C1).

Each test invokes the script as a subprocess so we test the real CLI
interface (exit codes, stderr output, JSON output).
"""

import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent.parent / 'bin' / \
    'parse_samplesheet.py'


def run(samplesheet, data_dir, tmp_path):
    out = tmp_path / 'out.json'
    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         '--samplesheet', str(samplesheet),
         '--data-dir', str(data_dir),
         '--out', str(out)],
        capture_output=True,
        text=True,
    )
    return result, out


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
    def test_valid_multirow_sheet_emits_in_order(self, tmp_path):
        d = make_data_dir(tmp_path, ['plant.fastq.gz', 'animal.fq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'PLANT-01,plant,plant.fastq.gz\n'
            'ANIMAL-01,animal,animal.fq.gz\n'
        ))
        result, out = run(csv, d, tmp_path)
        assert result.returncode == 0, result.stderr
        rows = json.loads(out.read_text())
        assert len(rows) == 2
        assert rows[0]['sample_id'] == 'PLANT-01'
        assert rows[1]['sample_id'] == 'ANIMAL-01'

    def test_valid_multifile_row_resolves_all_reads(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz', 'b.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'MULTI-01,animal,a.fastq.gz|b.fastq.gz\n'
        ))
        result, out = run(csv, d, tmp_path)
        assert result.returncode == 0, result.stderr
        rows = json.loads(out.read_text())
        assert len(rows[0]['reads']) == 2

    def test_reads_resolved_to_absolute_paths(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'S1,plant,a.fastq.gz\n'
        ))
        result, out = run(csv, d, tmp_path)
        assert result.returncode == 0, result.stderr
        rows = json.loads(out.read_text())
        assert Path(rows[0]['reads'][0]).is_absolute()

    def test_optional_columns_absent_from_header_ok(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'S1,plant,a.fastq.gz\n'
        ))
        result, out = run(csv, d, tmp_path)
        assert result.returncode == 0, result.stderr
        row = json.loads(out.read_text())[0]
        assert row['sample_info'] == ''
        assert row['sample_type'] == ''
        assert row['storage_location'] == ''

    def test_optional_columns_empty_cells_ok(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads,sample_info,sample_type,'
            'sample_receipt_date,storage_location\n'
            'S1,plant,a.fastq.gz,,,,\n'
        ))
        result, out = run(csv, d, tmp_path)
        assert result.returncode == 0, result.stderr

    def test_kingdom_case_normalised_to_lowercase(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'S1,Plant,a.fastq.gz\n'
        ))
        result, out = run(csv, d, tmp_path)
        assert result.returncode == 0, result.stderr
        assert json.loads(out.read_text())[0]['kingdom'] == 'plant'

    def test_bad_date_warns_but_passes_through(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads,sample_receipt_date\n'
            'S1,plant,a.fastq.gz,24/07/2026\n'
        ))
        result, out = run(csv, d, tmp_path)
        assert result.returncode == 0, result.stderr
        assert 'WARNING' in result.stderr
        assert json.loads(out.read_text())[0]['sample_receipt_date'] == \
            '24/07/2026'

    def test_all_extensions_accepted(self, tmp_path):
        exts = ['a.fastq', 'b.fastq.gz', 'c.fq', 'd.fq.gz']
        d = make_data_dir(tmp_path, exts)
        reads = '|'.join(exts)
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            f'S1,plant,{reads}\n'
        ))
        result, out = run(csv, d, tmp_path)
        assert result.returncode == 0, result.stderr
        assert len(json.loads(out.read_text())[0]['reads']) == 4


# ---------------------------------------------------------------------------
# Failure tests
# ---------------------------------------------------------------------------

class TestFailures:
    def test_duplicate_sample_id_fails(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'DUP-01,plant,a.fastq.gz\n'
            'DUP-01,animal,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path)
        assert result.returncode != 0
        assert 'DUP-01' in result.stderr

    def test_absolute_reads_path_fails(self, tmp_path):
        d = make_data_dir(tmp_path, [])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'S1,plant,/absolute/path.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path)
        assert result.returncode != 0
        assert 'absolute path' in result.stderr

    def test_missing_reads_file_fails(self, tmp_path):
        d = make_data_dir(tmp_path, [])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'S1,plant,nonexistent.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path)
        assert result.returncode != 0
        assert 'not found' in result.stderr

    def test_bad_kingdom_fails(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'S1,fungi,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path)
        assert result.returncode != 0
        assert 'fungi' in result.stderr

    def test_unknown_column_in_header_fails(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads,extra_column\n'
            'S1,plant,a.fastq.gz,value\n'
        ))
        result, _ = run(csv, d, tmp_path)
        assert result.returncode != 0
        assert 'extra_column' in result.stderr

    def test_missing_required_column_fails(self, tmp_path):
        d = make_data_dir(tmp_path, [])
        csv = write_csv(tmp_path, (
            'sample_id,reads\n'
            'S1,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path)
        assert result.returncode != 0
        assert 'kingdom' in result.stderr

    def test_bad_sample_id_pattern_fails(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'bad id!,plant,a.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path)
        assert result.returncode != 0

    def test_unsupported_reads_extension_fails(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.bam'])
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'S1,plant,a.bam\n'
        ))
        result, _ = run(csv, d, tmp_path)
        assert result.returncode != 0
        assert 'extension' in result.stderr

    def test_all_errors_collected_before_exit(self, tmp_path):
        d = make_data_dir(tmp_path, ['a.fastq.gz'])
        # 3 distinct problems: bad kingdom, duplicate id, absolute path
        csv = write_csv(tmp_path, (
            'sample_id,kingdom,reads\n'
            'S1,fungi,a.fastq.gz\n'
            'S1,plant,/abs/path.fastq.gz\n'
        ))
        result, _ = run(csv, d, tmp_path)
        assert result.returncode != 0
        assert result.stderr.count('ERROR:') >= 3
