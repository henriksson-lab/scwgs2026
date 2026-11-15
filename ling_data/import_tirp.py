#!/usr/bin/env python3
"""Import Ling's matched genomic/barcode/genomic FASTQ triplets as sorted BBGZ TIRP."""
import argparse
import collections
import csv
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import tempfile
import zlib

RAW_LIMIT = 60 * 1024
EOF = bytes([31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0,
             27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0])
ROLES = {"r1", "barcode", "r2_fragment"}


class DigestWriter:
    def __init__(self, handle):
        self.handle = handle
        self.sha = hashlib.sha256()
        self.size = 0

    def write(self, data):
        self.handle.write(data)
        self.sha.update(data)
        self.size += len(data)


class BBGZWriter:
    """One cell per block, BC last, mergeable sync-flushed raw DEFLATE."""
    def __init__(self, handle):
        self.out = DigestWriter(handle)
        self.cell = b""
        self.raw = bytearray()

    def record(self, cell, record):
        if not cell or len(cell) > 1024 or len(record) > RAW_LIMIT:
            raise ValueError("Cell/record exceeds BBGZ limits")
        if cell != self.cell or len(self.raw) + len(record) > RAW_LIMIT:
            self.flush()
            self.cell = cell
        self.raw.extend(record)

    def flush(self):
        if not self.raw:
            return
        compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
        compressed = compressor.compress(self.raw) + compressor.flush(zlib.Z_SYNC_FLUSH) + b"\x03\x00"
        xlen = 4 + len(self.cell) + 6
        size = 12 + xlen + len(compressed) + 8
        if size > 65536:
            raise ValueError("Compressed BBGZ block exceeds 64 KiB")
        self.out.write(b"\x1f\x8b\x08\x05" + struct.pack("<I", 0) + b"\x00\xff")
        self.out.write(struct.pack("<H", xlen))
        self.out.write(b"ID" + struct.pack("<H", len(self.cell)) + self.cell)
        self.out.write(b"BC\x02\x00" + struct.pack("<H", size - 1))
        self.out.write(compressed)
        self.out.write(struct.pack("<II", zlib.crc32(self.raw) & 0xffffffff, len(self.raw)))
        self.raw.clear()

    def finish(self):
        self.flush()
        self.out.write(EOF)


def fastq(path):
    with gzip.open(path, "rb") as handle:
        number = 0
        while header := handle.readline():
            number += 1
            seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
            if not seq or not plus or not qual or not header.startswith(b"@") or not plus.startswith(b"+"):
                raise ValueError(f"Malformed/truncated FASTQ {path}, record {number}")
            seq, qual = seq.rstrip(b"\r\n"), qual.rstrip(b"\r\n")
            if len(seq) != len(qual) or any(c in seq + qual for c in [9, 10, 13]):
                raise ValueError(f"Invalid FASTQ sequence/qualities {path}, record {number}")
            yield header.rstrip(b"\r\n"), seq, qual


def read_key(header):
    token = header[1:].split()[0].split(b"rbc:", 1)[0]
    # Deposited headers collapse the Illumina whitespace separating read
    # number from the final coordinate: ...:10491:N:0:0 vs ...:10492:N:0:0.
    match = re.fullmatch(rb"(.*)([12]):[YN]:[0-9]+:[^:]*", token)
    if match:
        return match[1]
    return token.removesuffix(b"/1").removesuffix(b"/2")


def source_identity(rows):
    return {role: {"filename": row["filename"], "md5": row["md5"], "bytes": int(row["bytes"])}
            for role, row in rows.items()}


def source_path(root, row):
    return root / "source" / row["library"] / row["filename"]


def verify_sources(root, groups):
    for _, rows in groups:
        for row in rows.values():
            path = source_path(root, row)
            receipt = path.with_name(path.name + ".verified.json")
            if not path.exists() or not receipt.exists():
                raise FileNotFoundError(f"Source download is not verified yet: {path}")
            saved = json.loads(receipt.read_text())
            if path.stat().st_size != int(row["bytes"]) or saved.get("md5") != row["md5"] \
                    or saved.get("bytes") != int(row["bytes"]):
                raise ValueError(f"Source/verification receipt disagrees with manifest: {path}")


def import_group(root, index, group, rows, sort_memory, sort_temp):
    library, barcode_a = group
    if not re.fullmatch(r"[a-z_]+", library) or not re.fullmatch(r"[ACGT]{8}", barcode_a):
        raise ValueError(f"Unsafe library/barcode-A identity: {group}")
    output = root / f"debarcoded.{index}.tirp.gz"
    histogram = output.with_name(output.name + ".hist")
    receipt = output.with_name(output.name + ".json")
    identity = {"library": library, "barcode_a": barcode_a, "sources": source_identity(rows)}
    if all(path.exists() for path in [output, histogram, receipt]):
        saved = json.loads(receipt.read_text())
        if all(saved.get(key) == value for key, value in identity.items()) \
                and saved.get("tirp_bytes") == output.stat().st_size \
                and saved.get("hist_bytes") == histogram.stat().st_size:
            print(f"Reusing {output.name}", flush=True)
            return saved
        raise RuntimeError(f"Existing output/receipt mismatch: {output}")
    if any(path.exists() for path in [output, histogram, receipt]):
        raise RuntimeError(f"Refusing incomplete existing output: {output}")
    group_temp = sort_temp / f"{library}_{barcode_a}"
    group_temp.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, LC_ALL="C")
    process = subprocess.Popen([
        "sort", "--stable", "--field-separator=\t", "--key=1,1", "--parallel=1",
        "--buffer-size=" + sort_memory, "--temporary-directory=" + str(group_temp),
    ], stdin=subprocess.PIPE, stdout=subprocess.PIPE, env=env)
    names, published = [], []
    try:
        iterators = {role: fastq(source_path(root, rows[role])) for role in ROLES}
        input_records = 0
        while True:
            records = {role: next(iterator, None) for role, iterator in iterators.items()}
            if all(record is None for record in records.values()):
                break
            if any(record is None for record in records.values()):
                raise ValueError(f"FASTQ triplet record counts differ: {group}")
            r1, bc, r2 = records["r1"], records["barcode"], records["r2_fragment"]
            if not read_key(r1[0]) == read_key(bc[0]) == read_key(r2[0]):
                raise ValueError(f"FASTQ triplet read IDs differ: {group}, record {input_records + 1}")
            if len(bc[1]) != 24 or not re.fullmatch(rb"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]{24}", bc[1]):
                raise ValueError(f"Expected 24-base observed barcode: {group}, record {input_records + 1}")
            cell = library.encode() + b"_" + barcode_a.encode() + b"_" + bc[1]
            # Seven separators preserve all eight fields including empty UMI.
            line = b"\t".join([cell, b"1", b"1", r1[1], r2[1], r1[2], r2[2], b""]) + b"\n"
            process.stdin.write(line)
            input_records += 1
        process.stdin.close()
        process.stdin = None
        data_fd, data_name = tempfile.mkstemp(prefix=".ling-tirp-", dir=root)
        names.append(Path(data_name))
        hist_fd, hist_name = tempfile.mkstemp(prefix=".ling-hist-", dir=root)
        names.append(Path(hist_name))
        with os.fdopen(data_fd, "wb") as data, os.fdopen(hist_fd, "wb") as hist:
            writer, hist_writer = BBGZWriter(data), DigestWriter(hist)
            current, count, cells, output_records = None, 0, 0, 0
            for line in process.stdout:
                fields = line.rstrip(b"\n").split(b"\t")
                if len(fields) != 8 or not line.endswith(b"\n"):
                    raise ValueError("Malformed sorted TIRP record")
                cell = fields[0]
                if current is not None and cell < current:
                    raise ValueError("External sort produced non-monotone cell IDs")
                if current != cell:
                    if current is not None:
                        hist_writer.write(current + b"\t" + str(count).encode() + b"\n")
                    current, count = cell, 0
                    cells += 1
                writer.record(cell, line)
                count += 1
                output_records += 1
            if current is not None:
                hist_writer.write(current + b"\t" + str(count).encode() + b"\n")
            writer.finish()
            for handle in [data, hist]:
                handle.flush()
                os.fsync(handle.fileno())
        status = process.wait()
        if status or input_records != output_records or input_records == 0:
            raise RuntimeError(f"External sort/import failed: status={status}, in={input_records}, out={output_records}")
        identity.update(records=input_records, cells=cells, tirp_bytes=writer.out.size,
                        hist_bytes=hist_writer.size, tirp_sha256=writer.out.sha.hexdigest(),
                        hist_sha256=hist_writer.sha.hexdigest())
        os.link(data_name, output)
        published.append(output)
        os.link(hist_name, histogram)
        published.append(histogram)
        # Receipt is written only after both complete data files are published.
        receipt_fd, receipt_name = tempfile.mkstemp(prefix=".ling-receipt-", dir=root)
        names.append(Path(receipt_name))
        with os.fdopen(receipt_fd, "w") as handle:
            json.dump(identity, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.link(receipt_name, receipt)
        published.append(receipt)
        print(f"Ready {output.name}: {input_records} read pairs, {cells} observed cells", flush=True)
        return identity
    except BaseException:
        if process.poll() is None:
            process.terminate()
        process.wait()
        for path in published:
            path.unlink(missing_ok=True)
        raise
    finally:
        for path in names:
            path.unlink(missing_ok=True)
        for iterator in locals().get("iterators", {}).values():
            iterator.close()
        if process.stdout:
            process.stdout.close()
        if process.stdin:
            try:
                process.stdin.close()
            except BrokenPipeError:
                pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/husky/henriksson/atrandi/v9_ling"))
    parser.add_argument("--manifest", type=Path, default=Path(__file__).with_name("manifest.tsv"))
    parser.add_argument("--sort-memory", default="4G", help="GNU sort buffer size; default 4G")
    parser.add_argument("--sort-temp", type=Path, help="GNU sort temporary storage; default ROOT/temp/ling_import")
    args = parser.parse_args()
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    sort_temp = args.sort_temp or root / "temp" / "ling_import"
    with (root / ".import.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        grouped = collections.defaultdict(dict)
        with args.manifest.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = row["library"], row["barcode_a"]
                if row["role"] in grouped[key]:
                    raise ValueError(f"Duplicate source role {key}, {row['role']}")
                grouped[key][row["role"]] = row
        groups = sorted(grouped.items())
        if not groups or any(set(rows) != ROLES for _, rows in groups):
            raise ValueError("Manifest must contain complete R1/barcode/R2-fragment triplets")
        verify_sources(root, groups)
        metadata = []
        for index, (group, rows) in enumerate(groups, 1):
            imported = import_group(root, index, group, rows, args.sort_memory, sort_temp)
            metadata.append({
                "prefix": group[0], "r1": rows["r1"]["filename"],
                "r2": rows["r2_fragment"]["filename"],
                "dir": str(source_path(root, rows["r1"]).parent), "barcode_a": group[1],
                "barcode_source": rows["barcode"]["filename"], "shard": index,
                "records": imported["records"], "cells": imported["cells"],
            })
        with tempfile.NamedTemporaryFile("w", dir=root, delete=False) as handle:
            temp = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=list(metadata[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(metadata)
        output = root / "debarcoded.meta"
        if output.exists():
            if output.read_bytes() != temp.read_bytes():
                raise RuntimeError(f"Refusing to overwrite different metadata: {output}")
            temp.unlink()
        else:
            os.link(temp, output)
            temp.unlink()
        print(f"Ready {len(groups)} sorted TIRPs and debarcoded.meta", flush=True)


if __name__ == "__main__":
    main()
