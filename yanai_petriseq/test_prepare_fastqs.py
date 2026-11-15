import gzip
import hashlib
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest

from prepare_fastqs import prepare_mate


@unittest.skipUnless(shutil.which("bgzip"), "bgzip is required")
class PrepareFastqsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.args = SimpleNamespace(raw_dir=root / "raw", output_dir=root / "bgzf",
                                    download=False, threads=1)
        self.args.raw_dir.mkdir()
        self.args.output_dir.mkdir()
        self.sequence = b"@read1\nACGTN\n+\nABCDE\n"
        self.source = self.args.raw_dir / "SRR123_1.fastq.gz"
        self.source.write_bytes(gzip.compress(self.sequence))
        self.row = dict(run_accession="SRR123", library_prefix="D1_lib1",
                        fastq_ftp="ftp.sra.ebi.ac.uk/x/SRR123_1.fastq.gz;ftp.sra.ebi.ac.uk/x/SRR123_2.fastq.gz",
                        fastq_md5=hashlib.md5(self.source.read_bytes()).hexdigest() + ";unused",
                        fastq_bytes=str(self.source.stat().st_size) + ";0")

    def test_recompression_preserves_every_byte_and_can_resume(self):
        output, reused = prepare_mate(self.row, 1, self.args)
        self.assertFalse(reused)
        self.assertEqual(output.name, "D1_lib1_S1_L001_R1_001.fastq.gz")
        self.assertEqual(gzip.decompress(output.read_bytes()), self.sequence)
        self.assertEqual(output.read_bytes()[:4], bytes([31, 139, 8, 4]))
        before = output.stat().st_mtime_ns
        again, reused = prepare_mate(self.row, 1, self.args)
        self.assertEqual(again, output)
        self.assertTrue(reused)
        self.assertEqual(before, output.stat().st_mtime_ns)

    def test_bad_source_checksum_is_rejected(self):
        self.row["fastq_md5"] = "0" * 32 + ";unused"
        with self.assertRaisesRegex(RuntimeError, "MD5 mismatch"):
            prepare_mate(self.row, 1, self.args)
        self.assertEqual(list(self.args.output_dir.iterdir()), [])

    def test_missing_receipt_never_silently_overwrites_output(self):
        output = self.args.output_dir / "D1_lib1_S1_L001_R1_001.fastq.gz"
        output.write_bytes(b"existing")
        with self.assertRaisesRegex(RuntimeError, "refusing overwrite"):
            prepare_mate(self.row, 1, self.args)
        self.assertEqual(output.read_bytes(), b"existing")


if __name__ == "__main__":
    unittest.main()
