#!/usr/bin/env python3
"""Compare an entire PETRI-seq library with an independent whitelist decoder.

Checks every input pair, every output record, lexical cell sorting, histogram
counts, and an order-independent SHA-256 multiset checksum of complete TIRP
records (including transcript quality and UMI). Memory scales with cell count.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import time

from verify_fastqs import barcode_tables


def record(stream):
    name = stream.readline()
    if not name:
        return None
    sequence = stream.readline().rstrip(b"\r\n")
    plus = stream.readline()
    quality = stream.readline().rstrip(b"\r\n")
    if not name.startswith(b"@") or not plus.startswith(b"+"):
        raise ValueError("Malformed FASTQ record")
    if not sequence or len(sequence) != len(quality):
        raise ValueError("Invalid sequence/quality length")
    return name.split()[0], sequence, quality


def decode(sequence, tables, prefix):
    sequence = sequence.upper()
    if len(sequence) < 57:
        return None
    for start, end, linker in ((14,29,b"GGTCCTTGGCTTCGC"),(36,50,b"CCTCCTACGCCAGA")):
        observed = sequence[start:end]
        if observed != linker and sum(a != b for a,b in zip(observed,linker)) > 2:
            return None
    a = tables[0].get(sequence[50:57])
    b = tables[1].get(sequence[29:36])
    c = tables[2].get(sequence[7:14])
    if a is None or b is None or c is None:
        return None
    return b"_".join((prefix,a,b,c))


def digest_value(line):
    return int.from_bytes(hashlib.sha256(line).digest(), "big")


def verify(args):
    started = time.monotonic()
    tables = [{k.encode(): (v[1].encode() if v[1] is not None else None)
               for k,v in table.items()} for table in barcode_tables(args.barcodes)]
    expected_cells = Counter()
    expected_digest = 0
    total = accepted = 0
    with gzip.open(args.r1,"rb") as left, gzip.open(args.r2,"rb") as right:
        while True:
            r1, r2 = record(left), record(right)
            if r1 is None and r2 is None:
                break
            if r1 is None or r2 is None or r1[0] != r2[0]:
                raise ValueError(f"Mates disagree after {total} pairs")
            total += 1
            cell = decode(r1[1],tables,args.library.encode())
            if cell is not None:
                accepted += 1
                expected_cells[cell] += 1
                line = b"\t".join((cell,b"1",b"1",b"",r2[1],b"",r2[2],r1[1][:7]))
                expected_digest += digest_value(line)
            if total % 1000000 == 0:
                print(f"Oracle: {total:,} pairs, {accepted:,} accepted",flush=True)
    if not accepted:
        raise ValueError("No accepted input pairs")

    observed_cells = Counter()
    observed_digest = 0
    observed = 0
    previous = b""
    with gzip.open(args.tirp,"rb") as stream:
        for raw in stream:
            line = raw.rstrip(b"\r\n")
            fields = line.split(b"\t")
            if len(fields)!=8:
                raise ValueError(f"TIRP record {observed+1} does not have eight columns")
            cell = fields[0]
            if cell < previous:
                raise ValueError(f"TIRP is not cell-sorted at record {observed+1}")
            previous = cell
            observed_cells[cell] += 1
            observed_digest += digest_value(line)
            observed += 1
            if observed % 1000000 == 0:
                print(f"TIRP: {observed:,} records checked",flush=True)
    modulus = 1 << 256
    if observed != accepted or observed_digest % modulus != expected_digest % modulus:
        raise ValueError("Output read multiset differs from independently decoded input (including cells, UMI, transcript and quality)")
    if observed_cells != expected_cells:
        raise ValueError("Output cell counts differ from independently decoded input")
    histogram = Counter()
    with args.hist.open("rb") as stream:
        for line in stream:
            cell,count = line.rstrip(b"\r\n").split(b"\t")
            if cell in histogram:
                raise ValueError("Duplicate histogram cell")
            histogram[cell] = int(count)
    if histogram != observed_cells:
        raise ValueError("Histogram differs from actual TIRP records")
    result = dict(status="passed",library=args.library,input_pairs=total,
                  accepted_pairs=accepted,rejected_pairs=total-accepted,
                  output_records=observed,barcode_combinations=len(observed_cells),
                  sorted_by_cell=True,all_records_checked=True,
                  histogram_matches=True,
                  multiset_sha256_sum=f"{expected_digest % modulus:064x}",
                  multiset_contents="All eight TIRP columns; original R2 sequence/quality and seven-base R1 UMI",
                  seconds=round(time.monotonic()-started,3),
                  r1=str(args.r1),r2=str(args.r2),tirp=str(args.tirp))
    args.report.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2),flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("r1","r2","tirp","hist","barcodes","report"):
        parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--library",required=True)
    verify(parser.parse_args())


if __name__ == "__main__":
    main()
