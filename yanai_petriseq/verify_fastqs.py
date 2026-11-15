#!/usr/bin/env python3
"""Audit prepared Yanai mate pairs and sample their PETRI-seq barcode layout."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import subprocess


def barcode_tables(path):
    """Independent nearest-whitelist oracle for up to one substitution."""
    rounds = [{}, {}, {}]
    with path.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            table = rounds[int(row["pos"]) - 1]
            sequence = row["seq"]
            for position in range(7):
                for base in "ACGTN":
                    variant = sequence[:position] + base + sequence[position + 1:]
                    distance = int(base != sequence[position])
                    candidate = (distance, row["well"])
                    previous = table.get(variant)
                    if previous is None or distance < previous[0]:
                        table[variant] = candidate
                    elif distance == previous[0] and candidate != previous:
                        table[variant] = (distance, None)
    return rounds


def barcode_match(sequence, rounds):
    if len(sequence) < 57:
        return False
    for observed, expected in [(sequence[14:29], "GGTCCTTGGCTTCGC"),
                               (sequence[36:50], "CCTCCTACGCCAGA")]:
        if sum(a != b for a, b in zip(observed, expected)) > 2:
            return False
    for (start, end), table in zip([(50, 57), (29, 36), (7, 14)], rounds):
        hit = table.get(sequence[start:end])
        if hit is None or hit[1] is None:
            return False
    return True


def read_record(stream):
    name = stream.readline()
    if not name:
        return None
    sequence, plus, quality = [stream.readline().rstrip("\r\n") for _ in range(3)]
    if not name.startswith("@") or not plus.startswith("+"):
        raise ValueError("Malformed FASTQ record")
    if len(sequence) != len(quality) or not sequence:
        raise ValueError("Invalid FASTQ sequence/quality length")
    return name.split()[0], sequence, quality


def verify_library(row, root, rounds, samples=1000, check_crc=False):
    result = dict(library_prefix=row["library_prefix"], run_accession=row["run_accession"],
                  status="waiting", records_checked=0, barcode_matches=0,
                  r1_lengths="", r2_lengths="", detail="")
    checksums = row["fastq_md5"].split(";")
    files = [root / "raw_bgzf" / f"{row['library_prefix']}_S1_L001_R{i}_001.fastq.gz"
             for i in [1, 2]]
    if not all(f.is_file() and f.with_name(f.name + ".json").is_file() for f in files):
        result["detail"] = "Waiting for both published FASTQs and verification receipts"
        return result
    try:
        for mate, path in enumerate(files, 1):
            receipt = json.loads(path.with_name(path.name + ".json").read_text())
            expected = dict(run=row["run_accession"], library=row["library_prefix"],
                            mate=mate, source_md5=checksums[mate - 1],
                            output_bytes=path.stat().st_size)
            if not all(receipt.get(key) == value for key, value in expected.items()):
                raise ValueError(f"Receipt mismatch: {path.name}")
            with path.open("rb") as stream:
                header = stream.read(18)
            if header[:4] != bytes([31, 139, 8, 4]) or header[12:16] != b"BC\x02\x00":
                raise ValueError(f"Expected standard BGZF header: {path.name}")
            if check_crc:
                subprocess.run(["bgzip", "-t", str(path)], check=True, capture_output=True)
        lengths = [set(), set()]
        with gzip.open(files[0], "rt", encoding="ascii") as left, \
                gzip.open(files[1], "rt", encoding="ascii") as right:
            for _ in range(samples):
                r1, r2 = read_record(left), read_record(right)
                if r1 is None and r2 is None:
                    break
                if r1 is None or r2 is None or r1[0] != r2[0]:
                    raise ValueError("R1/R2 read names or sampled record counts differ")
                result["records_checked"] += 1
                lengths[0].add(len(r1[1]))
                lengths[1].add(len(r2[1]))
                result["barcode_matches"] += int(barcode_match(r1[1].upper(), rounds))
        if result["records_checked"] == 0:
            raise ValueError("Empty FASTQs")
        if result["barcode_matches"] == 0:
            raise ValueError("No PETRI-seq barcodes match the Bascet layout in sampled R1")
        result.update(status="ready", r1_lengths=",".join(map(str, sorted(lengths[0]))),
                      r2_lengths=",".join(map(str, sorted(lengths[1]))),
                      detail="Receipts/BGZF headers/pair names/PETRI-seq sample validated"
                      + ("; full BGZF CRC passed" if check_crc else ""))
    except (ValueError, OSError, EOFError, subprocess.CalledProcessError) as error:
        result.update(status="invalid", detail=str(error))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--barcodes", type=Path,
                        default=Path("/home/mahogny/github/bascet/crates/bascet-cli/src/barcode/petriseq_barcodes.tsv"))
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--check-crc", action="store_true", help="Read and CRC-check every complete BGZF")
    parser.add_argument("--allow-incomplete", action="store_true", help="Permit pending downloads, never invalid files")
    parser.add_argument("--report", type=Path, help="Write per-library TSV report")
    args = parser.parse_args()
    if args.samples < 1:
        parser.error("--samples must be positive")
    root = args.root.resolve()
    with (root / "fastq_manifest.tsv").open() as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    if not manifest or len({row["library_prefix"] for row in manifest}) != len(manifest):
        parser.error("Manifest must contain unique, nonempty library prefixes")
    rounds = barcode_tables(args.barcodes)
    results = [verify_library(row, root, rounds, args.samples, args.check_crc) for row in manifest]
    if args.report:
        with args.report.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(results[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(results)
    counts = {status: sum(r["status"] == status for r in results)
              for status in ["ready", "waiting", "invalid"]}
    print(f"Libraries: {counts['ready']}/{len(results)} ready, "
          f"{counts['waiting']} waiting, {counts['invalid']} invalid")
    for result in results:
        if result["status"] == "invalid":
            print(f"{result['library_prefix']}: {result['detail']}")
    return int(counts["invalid"] > 0 or (counts["waiting"] > 0 and not args.allow_incomplete))


if __name__ == "__main__":
    raise SystemExit(main())
