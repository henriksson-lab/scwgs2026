import collections
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import struct
import tempfile
import unittest
import zlib

spec = importlib.util.spec_from_file_location("ling_import", Path(__file__).with_name("import_tirp.py"))
importer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(importer)


class ImportChecks(unittest.TestCase):
    def test_bbgz_blocks_and_mergeable_deflate_preserve_exact_bytes(self):
        out = io.BytesIO()
        writer = importer.BBGZWriter(out)
        cell = b"sewage_deep_ACTGACTG_" + b"ACGT" * 6
        line = b"\t".join([cell, b"1", b"1", b"A" * 140, b"C" * 115,
                           b"I" * 140, b"J" * 115, b""]) + b"\n"
        for _ in range(1000):
            writer.record(cell, line)
        writer.finish()
        packed = out.getvalue()
        self.assertEqual(gzip.decompress(packed), line * 1000)
        offset, payloads = 0, []
        while offset < len(packed):
            xlen = struct.unpack_from("<H", packed, offset + 10)[0]
            cursor, extras = offset + 12, {}
            while cursor < offset + 12 + xlen:
                name = packed[cursor:cursor + 2]
                size = struct.unpack_from("<H", packed, cursor + 2)[0]
                extras[name] = packed[cursor + 4:cursor + 4 + size]
                cursor += 4 + size
            block_size = struct.unpack("<H", extras[b"BC"])[0] + 1
            member = packed[offset:offset + block_size]
            raw = gzip.decompress(member)
            self.assertLessEqual(block_size, 65536)
            if raw:
                self.assertEqual(packed[offset + 3], 5)
                self.assertEqual(extras[b"ID"], cell)
                self.assertLessEqual(len(raw), 61440)
                self.assertTrue(raw.endswith(b"\n"))
                payloads.append(packed[cursor:offset + block_size - 8])
                self.assertTrue(payloads[-1].endswith(b"\x03\x00"))
            offset += block_size
        self.assertGreater(len(payloads), 1)
        merged = b"".join(payload[:-2] for payload in payloads) + b"\x03\x00"
        self.assertEqual(zlib.decompress(merged, -15), line * 1000)

    def fixture(self, root, mismatched=False):
        rows = {}
        barcodes = [b"N" * 24, b"ACGT" * 6, b"ACGT" * 6]
        records = {"r1": [], "barcode": [], "r2_fragment": []}
        for i, barcode in enumerate(barcodes):
            records["r1"].append((f"@read{i}/1".encode(), b"ACTG", b"ABCD"))
            records["r2_fragment"].append((f"@read{i}/2".encode(), b"GCA", b"EFG"))
            records["barcode"].append((f"@read{i}/2".encode(), barcode, b"I" * 24))
        if mismatched:
            records["barcode"][1] = (b"@different/2", b"ACGT" * 6, b"I" * 24)
        for role, reads in records.items():
            row = dict(library="sewage_deep", barcode_a="ACTGACTG", filename=role + ".fq.gz", role=role)
            path = importer.source_path(root, row)
            path.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(path, "wb") as stream:
                for header, sequence, quality in reads:
                    stream.write(header + b"\n" + sequence + b"\n+\n" + quality + b"\n")
            row.update(bytes=path.stat().st_size, md5=hashlib.md5(path.read_bytes()).hexdigest())
            path.with_name(path.name + ".verified.json").write_text(
                json.dumps({"bytes": row["bytes"], "md5": row["md5"]}))
            rows[role] = row
        return rows, barcodes

    def test_triplet_import_preserves_observed_barcode_and_measured_qualities(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows, barcodes = self.fixture(root)
            groups = [(("sewage_deep", "ACTGACTG"), rows)]
            importer.verify_sources(root, groups)
            receipt = importer.import_group(root, 1, groups[0][0], rows, "1M", root / "sort")
            self.assertEqual((receipt["records"], receipt["cells"]), (3, 2))
            lines = gzip.decompress((root / "debarcoded.1.tirp.gz").read_bytes()).splitlines()
            cells = [line.split(b"\t")[0] for line in lines]
            self.assertEqual(cells, sorted(cells))
            self.assertEqual(collections.Counter(cells),
                             collections.Counter(b"sewage_deep_ACTGACTG_" + bc for bc in barcodes))
            for line in lines:
                fields = line.split(b"\t")
                self.assertEqual(fields[1:], [b"1", b"1", b"ACTG", b"GCA", b"ABCD", b"EFG", b""])
            self.assertEqual(importer.import_group(root, 1, groups[0][0], rows, "1M", root / "sort"), receipt)

    def test_mismatched_triplet_ids_leave_no_published_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows, _ = self.fixture(root, mismatched=True)
            with self.assertRaisesRegex(ValueError, "read IDs differ"):
                importer.import_group(root, 1, ("sewage_deep", "ACTGACTG"), rows, "1M", root / "sort")
            self.assertFalse((root / "debarcoded.1.tirp.gz").exists())

    def test_ling_collapsed_header_mate_numbers_normalize_correctly(self):
        left = b"@NB502092:266:HHL37BGXT:1:11101:5649:10491:N:0:0rbc:ACGT"
        right = b"@NB502092:266:HHL37BGXT:1:11101:5649:10492:N:0:0rbc:ACGT"
        self.assertEqual(importer.read_key(left), importer.read_key(right))


if __name__ == "__main__":
    unittest.main()
