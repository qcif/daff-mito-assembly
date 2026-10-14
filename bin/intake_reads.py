#!/usr/bin/env python3
"""C1b: validate raw reads at intake and combine them into one FASTQ.

Checks every input file for gzip integrity and FASTQ record structure,
stopping each file at its first error. If every file passes, writes one
gzip FASTQ (gzip inputs byte-concatenated, plain inputs compressed).
Always writes an intake status file, and exits 0 once it has: a corrupt
input is data for the sample's report, not a Nextflow error.
"""

import argparse
import gzip
import json
import shutil
import sys
import zlib
from pathlib import Path

SCHEMA = 'wf5/intake/v1'
STATUS_VALID = 'valid'
STATUS_INVALID = 'invalid'

GZIP_MAGIC = b'\x1f\x8b'
COPY_CHUNK_BYTES = 1024 * 1024

ERR_GZIP_STREAM = 'gzip_stream'
ERR_TRUNCATED_RECORD = 'truncated_record'
ERR_LENGTH_MISMATCH = 'length_mismatch'
ERR_BAD_HEADER = 'bad_header'
ERR_BAD_SEPARATOR = 'bad_separator'
ERR_EMPTY = 'empty'

HEADER_PREFIX = b'@'
SEPARATOR_PREFIX = b'+'
LINES_PER_RECORD = 4
BYTES_PER_GB = 1_000_000_000
BYTES_PER_MB = 1_000_000

# Raised by the gzip/zlib layer on a truncated or corrupt stream.
GZIP_ERRORS = (EOFError, OSError, zlib.error)


def is_gzip(path: Path) -> bool:
    """True if the file starts with the gzip magic bytes."""
    with open(path, 'rb') as fh:
        return fh.read(len(GZIP_MAGIC)) == GZIP_MAGIC


def human_bytes(n: int) -> str:
    if n >= BYTES_PER_GB:
        return f'{n / BYTES_PER_GB:.2f} GB'
    return f'{n / BYTES_PER_MB:.1f} MB'


def make_error(kind: str, record: int, detail: str) -> dict:
    return {'kind': kind, 'record': record, 'detail': detail}


def scan_records(lines, state: dict):
    """Check FASTQ structure over an iterable of byte lines.

    `state` carries 'records' and 'bytes' counters, updated as records
    pass so they survive an exception raised by `lines`. Returns None,
    or an error dict at the first bad record without pulling further.
    """
    it = iter(lines)
    while True:
        header = next(it, None)
        if header is None:
            break
        record = state['records'] + 1
        rest = [next(it, None) for _ in range(LINES_PER_RECORD - 1)]
        if not header.startswith(HEADER_PREFIX):
            return make_error(
                ERR_BAD_HEADER, record,
                'header line does not start with "@"')
        if None in rest:
            return make_error(
                ERR_TRUNCATED_RECORD, record,
                'file ends partway through a record')
        seq, sep, qual = rest
        if not sep.startswith(SEPARATOR_PREFIX):
            return make_error(
                ERR_BAD_SEPARATOR, record,
                'separator line does not start with "+"')
        if len(seq.rstrip(b'\r\n')) != len(qual.rstrip(b'\r\n')):
            return make_error(
                ERR_LENGTH_MISMATCH, record,
                'sequence and quality lengths differ')
        state['records'] += 1
        state['bytes'] += (
            len(header) + len(seq) + len(sep) + len(qual))
    if state['records'] == 0:
        return make_error(ERR_EMPTY, 0, 'file contains no reads')
    return None


def check_file(path: Path) -> dict:
    """Validate one input file and return its intake record."""
    compressed = is_gzip(path)
    opener = gzip.open if compressed else open
    state = {'records': 0, 'bytes': 0}
    try:
        with opener(path, 'rb') as fh:
            error = scan_records(fh, state)
    except GZIP_ERRORS as exc:
        error = make_error(
            ERR_GZIP_STREAM, state['records'] + 1,
            f'invalid compressed data after '
            f'{human_bytes(state["bytes"])} decompressed ({exc})')
    return {
        'name': path.name,
        'compressed': compressed,
        'records': state['records'],
        'bytes': path.stat().st_size,
        'error': error,
    }


def write_reads(paths: list, out_path: Path) -> None:
    """Combine inputs into one multi-member gzip file."""
    with open(out_path, 'wb') as out:
        for path in paths:
            if is_gzip(path):
                with open(path, 'rb') as src:
                    shutil.copyfileobj(src, out, COPY_CHUNK_BYTES)
            else:
                with open(path, 'rb') as src:
                    with gzip.GzipFile(fileobj=out, mode='wb') as gz:
                        shutil.copyfileobj(src, gz, COPY_CHUNK_BYTES)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reads', nargs='+', required=True, type=Path)
    parser.add_argument('--out-reads', required=True, type=Path)
    parser.add_argument('--out-status', required=True, type=Path)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    files = [check_file(path) for path in args.reads]
    valid = all(f['error'] is None for f in files)
    status = {
        '$schema': SCHEMA,
        'status': STATUS_VALID if valid else STATUS_INVALID,
        'files': files,
        'total_records': sum(f['records'] for f in files),
    }
    if valid:
        write_reads(args.reads, args.out_reads)
    args.out_status.write_text(json.dumps(status, indent=2) + '\n')
    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
