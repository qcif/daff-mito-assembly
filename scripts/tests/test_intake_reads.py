"""Unit tests for bin/intake_reads.py (C1b) — task 61 §2.1.

Real gzip/plain files built in a temp dir; no mocks.
"""

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN_DIR))

import intake_reads  # noqa: E402


def make_records(n: int, prefix: str = "r") -> bytes:
    return b"".join(
        b"@%s%d\nACGTACGT\n+\nIIIIIIII\n" % (prefix.encode(), i)
        for i in range(n)
    )


class IntakeCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = Path(self._tmp.name)

    def write(self, name: str, data: bytes) -> Path:
        path = self.tmp / name
        path.write_bytes(data)
        return path

    def write_gz(self, name: str, data: bytes) -> Path:
        return self.write(name, gzip.compress(data))

    def run_main(self, *inputs: Path) -> dict:
        out = self.tmp / "out.fastq.gz"
        status = self.tmp / "intake.json"
        rc = intake_reads.main([
            "--reads", *map(str, inputs),
            "--out-reads", str(out),
            "--out-status", str(status),
        ])
        self.assertEqual(rc, 0)
        return json.loads(status.read_text())

    def first_error(self, *inputs: Path) -> dict:
        result = self.run_main(*inputs)
        self.assertEqual(result["status"], "invalid")
        return result["files"][-1]["error"]


class TestValidInputs(IntakeCase):
    def test_gzip_input_is_valid_and_copied_byte_identical(self):
        src = self.write_gz("a.fastq.gz", make_records(5))
        result = self.run_main(src)
        self.assertEqual(result["$schema"], "wf5/intake/v1")
        self.assertEqual(result["status"], "valid")
        self.assertEqual(result["total_records"], 5)
        entry = result["files"][0]
        self.assertEqual(entry["name"], "a.fastq.gz")
        self.assertTrue(entry["compressed"])
        self.assertEqual(entry["records"], 5)
        self.assertEqual(entry["bytes"], src.stat().st_size)
        self.assertIsNone(entry["error"])
        self.assertEqual(
            (self.tmp / "out.fastq.gz").read_bytes(), src.read_bytes())

    def test_plain_input_is_compressed(self):
        src = self.write("a.fastq", make_records(3))
        result = self.run_main(src)
        self.assertEqual(result["status"], "valid")
        self.assertFalse(result["files"][0]["compressed"])
        out = gzip.decompress((self.tmp / "out.fastq.gz").read_bytes())
        self.assertEqual(out, make_records(3))

    def test_multi_member_gzip_is_valid(self):
        data = (
            gzip.compress(make_records(2, "a"))
            + gzip.compress(make_records(3, "b")))
        src = self.write("m.fastq.gz", data)
        result = self.run_main(src)
        self.assertEqual(result["status"], "valid")
        self.assertEqual(result["files"][0]["records"], 5)

    def test_mixed_inputs_give_ordered_gzip_output(self):
        gz = self.write_gz("a.fastq.gz", make_records(2, "a"))
        plain = self.write("b.fastq", make_records(2, "b"))
        result = self.run_main(gz, plain)
        self.assertEqual(result["status"], "valid")
        self.assertEqual(result["total_records"], 4)
        out = gzip.decompress((self.tmp / "out.fastq.gz").read_bytes())
        self.assertEqual(
            out, make_records(2, "a") + make_records(2, "b"))

    def test_gzip_detected_by_magic_not_extension(self):
        gz_named_plain = self.write_gz("a.fastq", make_records(2))
        plain_named_gz = self.write("b.fastq.gz", make_records(2))
        result = self.run_main(gz_named_plain, plain_named_gz)
        self.assertEqual(result["status"], "valid")
        self.assertEqual(
            [f["compressed"] for f in result["files"]], [True, False])

    def test_crlf_line_endings_are_valid(self):
        src = self.write(
            "a.fastq", b"@r1\r\nACGT\r\n+\r\nIIII\r\n")
        self.assertEqual(self.run_main(src)["status"], "valid")

    def test_missing_final_newline_is_valid(self):
        src = self.write("a.fastq", b"@r1\nACGT\n+\nIIII")
        self.assertEqual(self.run_main(src)["status"], "valid")


class TestErrors(IntakeCase):
    def test_truncated_gzip(self):
        data = gzip.compress(make_records(20000))
        src = self.write("a.fastq.gz", data[:len(data) // 2])
        error = self.first_error(src)
        self.assertEqual(error["kind"], "gzip_stream")
        self.assertGreater(error["record"], 1)
        self.assertIn("invalid compressed data after", error["detail"])

    def test_corrupt_gzip(self):
        data = bytearray(gzip.compress(make_records(20000)))
        data[len(data) // 2] ^= 0xFF
        src = self.write("a.fastq.gz", bytes(data))
        # Highly repetitive test data may decode to garbage records
        # before the CRC check trips, so only assert it is rejected.
        self.assertIsNotNone(self.first_error(src))

    def test_trailing_garbage_after_gzip(self):
        data = gzip.compress(make_records(3)) + b"garbage-bytes"
        src = self.write("a.fastq.gz", data)
        result = self.run_main(src)
        entry = result["files"][0]
        self.assertEqual(result["status"], "invalid")
        self.assertEqual(entry["error"]["kind"], "gzip_stream")
        self.assertEqual(entry["error"]["record"], 4)
        self.assertEqual(entry["records"], 3)

    def test_gigabyte_sized_detail(self):
        self.assertEqual(
            intake_reads.human_bytes(1_980_000_000), "1.98 GB")
        self.assertEqual(intake_reads.human_bytes(2_500_000), "2.5 MB")

    def test_length_mismatch(self):
        data = make_records(2) + b"@bad\nACGT\n+\nIII\n"
        error = self.first_error(self.write("a.fastq", data))
        self.assertEqual(error["kind"], "length_mismatch")
        self.assertEqual(error["record"], 3)

    def test_bad_header(self):
        data = make_records(1) + b"r2\nACGT\n+\nIIII\n"
        error = self.first_error(self.write("a.fastq", data))
        self.assertEqual(error["kind"], "bad_header")
        self.assertEqual(error["record"], 2)

    def test_bad_separator(self):
        data = make_records(1) + b"@r2\nACGT\n-\nIIII\n"
        error = self.first_error(self.write("a.fastq", data))
        self.assertEqual(error["kind"], "bad_separator")
        self.assertEqual(error["record"], 2)

    def test_final_record_missing_quality(self):
        data = make_records(2) + b"@r3\nACGT\n+\n"
        error = self.first_error(self.write("a.fastq", data))
        self.assertEqual(error["kind"], "truncated_record")
        self.assertEqual(error["record"], 3)

    def test_empty_file(self):
        error = self.first_error(self.write("a.fastq", b""))
        self.assertEqual(error["kind"], "empty")

    def test_gzip_file_with_empty_content(self):
        error = self.first_error(self.write_gz("a.fastq.gz", b""))
        self.assertEqual(error["kind"], "empty")

    def test_second_file_bad_first_still_counted_no_output(self):
        good = self.write_gz("a.fastq.gz", make_records(4))
        bad = self.write("b.fastq", b"@r\nACGT\n+\nII\n")
        result = self.run_main(good, bad)
        self.assertEqual(result["status"], "invalid")
        self.assertEqual(result["files"][0]["records"], 4)
        self.assertIsNone(result["files"][0]["error"])
        self.assertEqual(
            result["files"][1]["error"]["kind"], "length_mismatch")
        self.assertEqual(result["total_records"], 4)
        self.assertFalse((self.tmp / "out.fastq.gz").exists())

    def test_several_bad_files_are_all_reported(self):
        a = self.write("a.fastq", b"")
        b = self.write("b.fastq", b"junk\n")
        result = self.run_main(a, b)
        self.assertEqual(
            [f["error"]["kind"] for f in result["files"]],
            ["empty", "bad_header"])


class TestScanRecords(unittest.TestCase):
    def test_stops_pulling_after_first_error(self):
        def lines():
            yield b"@r1\n"
            yield b"ACGT\n"
            yield b"+\n"
            yield b"II\n"  # length mismatch: stop here
            raise AssertionError("reader consumed past the first error")

        state = {"records": 0, "bytes": 0}
        error = intake_reads.scan_records(lines(), state)
        self.assertEqual(error["kind"], "length_mismatch")
        self.assertEqual(state["records"], 0)

    def test_state_survives_exception_from_source(self):
        def lines():
            yield from make_records(2).splitlines(keepends=True)
            raise EOFError("boom")

        state = {"records": 0, "bytes": 0}
        with self.assertRaises(EOFError):
            intake_reads.scan_records(lines(), state)
        self.assertEqual(state["records"], 2)


class TestCli(unittest.TestCase):
    def test_missing_arguments_exit_nonzero(self):
        with self.assertRaises(SystemExit) as ctx:
            with open("/dev/null", "w") as devnull:
                real_stderr, sys.stderr = sys.stderr, devnull
                try:
                    intake_reads.main([])
                finally:
                    sys.stderr = real_stderr
        self.assertNotEqual(ctx.exception.code, 0)


if __name__ == "__main__":
    unittest.main()
