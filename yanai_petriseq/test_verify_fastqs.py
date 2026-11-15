import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from verify_fastqs import barcode_tables, verify_library


@unittest.skipUnless(shutil.which("bgzip"), "bgzip is required")
class VerifyFastqsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "raw_bgzf").mkdir()
        barcodes = self.root / "barcodes.tsv"
        barcodes.write_text("pos\twell\tseq\n1\ta1\tCAGAGAA\n2\tb1\tGACCTTA\n3\tc1\tGTGTGAA\n")
        self.rounds = barcode_tables(barcodes)
        self.row = dict(library_prefix="D1_lib1", run_accession="SRR123", fastq_md5="abc;def")
        self.r1 = "ACGTNACGTGTGAAGGTCCTTGGCTTCGCGACCTTACCTCCTACGCCAGACAGAGAAT"
        self.files = []
        for mate, sequence in [(1, self.r1), (2, "ACGT" * 6 + "AC")]:
            path = self.root / "raw_bgzf" / f"D1_lib1_S1_L001_R{mate}_001.fastq.gz"
            record = f"@read1\n{sequence}\n+\n{'I' * len(sequence)}\n".encode()
            packed = subprocess.run(["bgzip", "-c"], input=record,
                                    capture_output=True, check=True).stdout
            path.write_bytes(packed)
            path.with_name(path.name + ".json").write_text(json.dumps(dict(
                run="SRR123", library="D1_lib1", mate=mate,
                source_md5=["abc", "def"][mate - 1], output_bytes=len(packed)
            )))
            self.files.append(path)

    def test_pair_names_layout_and_full_crc_are_validated(self):
        result = verify_library(self.row, self.root, self.rounds, check_crc=True)
        self.assertEqual(result["status"], "ready", result)
        self.assertEqual(result["records_checked"], 1)
        self.assertEqual(result["barcode_matches"], 1)
        self.assertEqual((result["r1_lengths"], result["r2_lengths"]), ("58", "26"))

    def test_mismatched_pair_names_fail(self):
        sequence = "ACGT" * 6 + "AC"
        packed = subprocess.run(["bgzip", "-c"],
            input=f"@different\n{sequence}\n+\n{'I' * len(sequence)}\n".encode(),
            capture_output=True, check=True).stdout
        self.files[1].write_bytes(packed)
        receipt_path = self.files[1].with_name(self.files[1].name + ".json")
        receipt = json.loads(receipt_path.read_text())
        receipt["output_bytes"] = len(packed)
        receipt_path.write_text(json.dumps(receipt))
        result = verify_library(self.row, self.root, self.rounds)
        self.assertEqual(result["status"], "invalid")
        self.assertIn("read names", result["detail"])

    def test_missing_mate_is_waiting_and_incorrect_receipt_is_invalid(self):
        receipt_path = self.files[1].with_name(self.files[1].name + ".json")
        receipt = receipt_path.read_text()
        receipt_path.unlink()
        self.assertEqual(verify_library(self.row, self.root, self.rounds)["status"], "waiting")
        incorrect = json.loads(receipt)
        incorrect["run"] = "SRR999"
        receipt_path.write_text(json.dumps(incorrect))
        self.assertEqual(verify_library(self.row, self.root, self.rounds)["status"], "invalid")


if __name__ == "__main__":
    unittest.main()
