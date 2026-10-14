"""Tests for the Cloudgene post-run scripts in deploy/cloudgene/."""

import json
import subprocess
import zipfile
from pathlib import Path

DEPLOY_DIR = Path(__file__).resolve().parents[2] / 'deploy' / 'cloudgene'
COLLECT_ERRORS = DEPLOY_DIR / 'collect_errors.sh'
ZIP_REPORTS = DEPLOY_DIR / 'zip_reports.sh'

STATUSES = ['ok', 'error', 'invalid_input']


def make_outdir(tmp_path, samples, manifest=True):
    """Build a fake outdir; samples maps sample_id -> (status, reason)."""
    if manifest:
        breakdown = {s: 0 for s in STATUSES}
        for status, _ in samples.values():
            breakdown[status] += 1
        (tmp_path / 'run_manifest.json').write_text(json.dumps({
            'status_breakdown': breakdown,
            'samples': [
                {
                    'sample_id': sample_id,
                    'sample_status': status,
                    'sample_status_reason': reason,
                }
                for sample_id, (status, reason) in samples.items()
            ],
        }))
    return tmp_path


def write_log(outdir, process, sample_id, ext, text):
    path = outdir / 'errors' / process / f'{sample_id}.{ext}'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def zip_files(path):
    with zipfile.ZipFile(path) as zf:
        return sorted(n for n in zf.namelist() if not n.endswith('/'))


def run(script, outdir):
    return subprocess.run(
        ['bash', str(script), '--dir', str(outdir)],
        capture_output=True, text=True)


def test_no_errors_dir(tmp_path):
    outdir = make_outdir(tmp_path, {'S1': ('ok', None)})
    result = run(COLLECT_ERRORS, outdir)
    assert result.returncode == 0
    assert 'No task errors' in result.stdout
    assert 'ok: 1' in result.stdout
    assert not (outdir / 'errors.zip').exists()


def test_error_logs_grouped_by_sample(tmp_path):
    outdir = make_outdir(tmp_path, {
        'S1': ('error', 'FILTLONG failed'),
        'S2': ('error', 'METAFLYE failed'),
    })
    write_log(outdir, 'FILTLONG', 'S1', 'err', 'filtlong boom\n')
    write_log(outdir, 'METAFLYE', 'S2', 'err', 'flye boom\n')
    write_log(outdir, 'METAFLYE', 'S2', 'out', 'ignored stdout\n')
    result = run(COLLECT_ERRORS, outdir)
    assert result.returncode == 0
    out = result.stdout
    assert 'S1 [error]: FILTLONG failed' in out
    assert out.index('=== S1 ===') < out.index('filtlong boom')
    assert out.index('=== S2 ===') < out.index('flye boom')
    assert out.index('filtlong boom') < out.index('=== S2 ===')
    assert 'ignored stdout' not in out
    assert zip_files(outdir / 'errors.zip') == [
        'errors/FILTLONG/S1.err',
        'errors/METAFLYE/S2.err',
        'errors/METAFLYE/S2.out',
    ]


def test_tail_limited_to_ten_lines(tmp_path):
    outdir = make_outdir(tmp_path, {'S1': ('error', 'x')})
    write_log(outdir, 'FILTLONG', 'S1', 'err',
              ''.join(f'line{i}\n' for i in range(30)))
    out = run(COLLECT_ERRORS, outdir).stdout
    assert 'line20' in out and 'line29' in out
    assert 'line19' not in out


def test_invalid_input_in_summary_without_logs(tmp_path):
    outdir = make_outdir(tmp_path, {
        'S1': ('ok', None),
        'S2': ('invalid_input', 'reads file is truncated'),
    })
    result = run(COLLECT_ERRORS, outdir)
    assert result.returncode == 0
    assert 'S2 [invalid_input]: reads file is truncated' in result.stdout
    assert 'invalid_input: 1' in result.stdout


def test_collect_errors_missing_manifest(tmp_path):
    result = run(COLLECT_ERRORS, make_outdir(tmp_path, {}, manifest=False))
    assert result.returncode != 0
    assert 'run_manifest.json not found' in result.stdout


def make_reports(outdir, with_report):
    (outdir / 'run-report.html').write_text('<html>run</html>')
    (outdir / 'run_manifest.json').write_text('{}')
    for sample_id in ('S1', 'S2', 'S3'):
        (outdir / sample_id).mkdir()
        (outdir / sample_id / 'metadata.json').write_text('{}')
        if sample_id in with_report:
            (outdir / sample_id / 'report.html').write_text(sample_id)


def test_zip_reports_layout(tmp_path):
    make_reports(tmp_path, with_report=['S1', 'S3'])
    result = run(ZIP_REPORTS, tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert zip_files(tmp_path / 'reports.zip') == [
        'reports/S1/report.html',
        'reports/S3/report.html',
        'reports/run-report.html',
        'reports/run_manifest.json',
    ]
    with zipfile.ZipFile(tmp_path / 'reports.zip') as zf:
        assert zf.read('reports/S3/report.html') == b'S3'


def test_zip_reports_idempotent(tmp_path):
    make_reports(tmp_path, with_report=['S1', 'S2'])
    assert run(ZIP_REPORTS, tmp_path).returncode == 0
    with zipfile.ZipFile(tmp_path / 'reports.zip') as zf:
        first = {n: zf.read(n) for n in zf.namelist()}
    assert run(ZIP_REPORTS, tmp_path).returncode == 0
    with zipfile.ZipFile(tmp_path / 'reports.zip') as zf:
        assert {n: zf.read(n) for n in zf.namelist()} == first
    assert not (tmp_path / 'reports').exists()
    assert sorted(p.name for p in tmp_path.iterdir() if p.is_file()) == [
        'reports.zip', 'run-report.html', 'run_manifest.json']


def test_zip_reports_missing_run_report(tmp_path):
    (tmp_path / 'S1').mkdir()
    (tmp_path / 'S1' / 'report.html').write_text('x')
    result = run(ZIP_REPORTS, tmp_path)
    assert result.returncode != 0
    assert not (tmp_path / 'reports.zip').exists()
