#!/usr/bin/env python3
"""Assertions for the stub-run failure-injection cases (task 61 §3.2).

Usage: check_injection.py CASE OUTDIR CONSOLE_LOG TASKS_TSV

CASE is one of: nanoplot, mixed, all_fail. TASKS_TSV is the output of
`nextflow log last -f name,status`.
"""

import json
import sys
from pathlib import Path

SAMPLES = ['STUB-01', 'STUB-02']
WARNING = 'WF5 WARNING:'
QC_UNAVAILABLE = 'NanoPlot read QC is unavailable'


def load_tasks(path):
    tasks = {}
    for line in Path(path).read_text().splitlines():
        name, _, status = line.rpartition('\t')
        tasks.setdefault(name.split(' (')[0], []).append(
            (name, status))
    return tasks


def main(case, outdir, console_log, tasks_tsv):
    out = Path(outdir)
    errors = []
    console = Path(console_log).read_text()
    tasks = load_tasks(tasks_tsv)

    if WARNING not in console:
        errors.append(f'Console output lacks the {WARNING!r} summary')

    manifest = json.loads((out / 'run_manifest.json').read_text())
    by_id = {s['sample_id']: s for s in manifest['samples']}
    if sorted(by_id) != SAMPLES:
        errors.append(f'run_manifest.json samples: {sorted(by_id)}')

    collate = [name for name, _ in tasks.get('COLLATE', [])]
    if case == 'all_fail':
        if collate:
            errors.append(f'COLLATE ran unexpectedly: {collate}')
        for sample_id in SAMPLES:
            if by_id.get(sample_id, {}).get('sample_status') != 'error':
                errors.append(f'{sample_id} is not status error')
        if manifest['status_breakdown']['error'] != len(SAMPLES):
            errors.append('status_breakdown.error != 2')
    else:
        if sorted(collate) != [f'COLLATE ({s})' for s in SAMPLES]:
            errors.append(f'COLLATE should run once per sample: {collate}')
        for sample_id in SAMPLES:
            for fname in ('metadata.json', 'report.html'):
                if not (out / sample_id / fname).exists():
                    errors.append(f'Missing {sample_id}/{fname}')
        meta1 = json.loads((out / 'STUB-01' / 'metadata.json').read_text())
        meta2 = json.loads((out / 'STUB-02' / 'metadata.json').read_text())
        if meta1['read_qc']['unavailable_reason'] is not None:
            errors.append('STUB-01 QC unexpectedly unavailable')
        if QC_UNAVAILABLE not in (meta2['read_qc']['unavailable_reason']):
            errors.append('STUB-02 QC section not marked unavailable')
        expected = 'fail' if case == 'mixed' else 'ok'
        if meta2['sample_status'] != expected:
            errors.append(f'STUB-02 sample_status != {expected}')

    if errors:
        for e in errors:
            print(f'ERROR: {e}', file=sys.stderr)
        return 1
    print(f'Failure injection {case} OK')
    return 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:5]))
