import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("ling_download", Path(__file__).with_name("download.py"))
download = importlib.util.module_from_spec(spec)
spec.loader.exec_module(download)


class DownloadChecks(unittest.TestCase):
    def row(self, content):
        return {
            "record": "13323680", "library": "sewage_deep", "barcode_a": "ATGTCTGC",
            "role": "barcode", "filename": "fixture.fastq.gz", "bytes": str(len(content)),
            "md5": hashlib.md5(content).hexdigest(), "url": "https://zenodo.org/fixture",
        }

    def test_existing_data_must_match_pinned_checksum(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            row = self.row(b"expected")
            source = root / "source" / row["library"] / row["filename"]
            source.parent.mkdir(parents=True)
            source.write_bytes(b"wrong!!!")
            with self.assertRaisesRegex(RuntimeError, "size/MD5"):
                download.fetch(row, root)
            self.assertFalse(source.with_name(source.name + ".verified.json").exists())
            self.assertEqual(source.read_bytes(), b"wrong!!!")

    def test_verified_existing_data_are_reused_with_receipt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            row = self.row(b"expected")
            source = root / "source" / row["library"] / row["filename"]
            source.parent.mkdir(parents=True)
            source.write_bytes(b"expected")
            self.assertEqual(download.fetch(row, root), source)
            receipt = json.loads(source.with_name(source.name + ".verified.json").read_text())
            self.assertEqual(receipt["md5"], row["md5"])
            self.assertEqual(receipt["bytes"], 8)

    def test_manifest_has_32_complete_triplets_and_96_unique_files(self):
        with Path(__file__).with_name("manifest.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 96)
        self.assertEqual(len({(r["record"], r["filename"]) for r in rows}), 96)
        groups = {}
        for row in rows:
            groups.setdefault((row["library"], row["barcode_a"]), set()).add(row["role"])
            self.assertTrue(row["url"].startswith("https://zenodo.org/api/records/"))
            self.assertEqual(len(row["md5"]), 32)
            self.assertGreater(int(row["bytes"]), 0)
        self.assertEqual(len(groups), 32)
        self.assertTrue(all(roles == {"r1", "barcode", "r2_fragment"} for roles in groups.values()))
        self.assertEqual(sum(int(r["bytes"]) for r in rows), 155063227653)


if __name__ == "__main__":
    unittest.main()
