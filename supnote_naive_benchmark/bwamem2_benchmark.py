#!/usr/bin/env python3
"""Generate a 10k-pair test set and benchmark bwa-mem2."""

import csv
import gzip
import shutil
import statistics
import subprocess
import time
from pathlib import Path


# Edit these settings if needed, then run: python3 bwamem2_benchmark.py
REFERENCE = Path("/husky/henriksson/atrandi/bwamem_ref/human_phi_xantham/all.fa")
SOURCE_R1 = Path(
    "/husky/henriksson/atrandi/rawdata/260116_mock_pta/filtered/"
    "fmock_pta2_R1_001.fastq.gz"
)
SOURCE_R2 = Path(
    "/husky/henriksson/atrandi/rawdata/260116_mock_pta/filtered/"
    "fmock_pta2_R2_001.fastq.gz"
)
BWA_MEM2 = "bwa-mem2"
PAIRS = 10_000
READ_LENGTH = 150
TOTAL_PAIRS = 800_000_000
THREADS = 10
REPEATS = 5
WARMUPS = 1

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
R1 = DATA_DIR / f"reads_{PAIRS}_R1.fastq.gz"
R2 = DATA_DIR / f"reads_{PAIRS}_R2.fastq.gz"


def read_record(handle, source):
    record = [handle.readline() for _ in range(4)]
    if any(line == "" for line in record):
        raise RuntimeError(f"{source} contains fewer than {PAIRS} records")
    return record


def read_name(header):
    return header.removeprefix("@").split()[0].removesuffix("/1").removesuffix("/2")


def generate_reads():
    """Copy the first PAIRS mates, trimming sequences and qualities to 150 bp."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Generating {PAIRS:,} pairs of {READ_LENGTH} bp reads...")
    with (
        gzip.open(SOURCE_R1, "rt") as source1,
        gzip.open(SOURCE_R2, "rt") as source2,
        gzip.open(R1, "wt", compresslevel=6) as output1,
        gzip.open(R2, "wt", compresslevel=6) as output2,
    ):
        for pair_number in range(1, PAIRS + 1):
            mate1 = read_record(source1, SOURCE_R1)
            mate2 = read_record(source2, SOURCE_R2)
            if read_name(mate1[0]) != read_name(mate2[0]):
                raise RuntimeError(f"mate names differ at pair {pair_number}")
            for mate, output in ((mate1, output1), (mate2, output2)):
                sequence = mate[1].rstrip("\r\n")
                quality = mate[3].rstrip("\r\n")
                if len(sequence) < READ_LENGTH or len(quality) < READ_LENGTH:
                    raise RuntimeError(f"pair {pair_number} is shorter than {READ_LENGTH} bp")
                output.writelines(
                    (mate[0], sequence[:READ_LENGTH] + "\n", mate[2], quality[:READ_LENGTH] + "\n")
                )
    print(f"Wrote {R1} and {R2}")


def validate_inputs():
    bwa = shutil.which(BWA_MEM2)
    if bwa is None:
        raise RuntimeError(f"cannot find {BWA_MEM2} on PATH")
    for path in (REFERENCE, SOURCE_R1, SOURCE_R2):
        if not path.is_file():
            raise RuntimeError(f"missing input: {path}")
    for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"):
        index = Path(f"{REFERENCE}{suffix}")
        if not index.is_file():
            raise RuntimeError(f"missing bwa-mem2 index: {index}")
    if TOTAL_PAIRS % PAIRS:
        raise RuntimeError("TOTAL_PAIRS must be divisible by PAIRS")
    return bwa


def timed_run(bwa, label, repeat, read1, read2):
    command = [
        bwa, "mem", "-t", str(THREADS), str(REFERENCE), str(read1), str(read2)
    ]
    log = RESULTS_DIR / "logs" / f"{label}_{repeat}.log"
    start = time.perf_counter()
    with log.open("w") as stderr:
        completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=stderr)
    elapsed = time.perf_counter() - start
    if completed.returncode:
        raise RuntimeError(f"{label} failed; see {log}")
    return elapsed


def benchmark(bwa):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "logs").mkdir(exist_ok=True)
    empty1 = RESULTS_DIR / "empty_R1.fastq"
    empty2 = RESULTS_DIR / "empty_R2.fastq"
    empty1.write_text("")
    empty2.write_text("")

    print(f"Warming the filesystem cache ({WARMUPS} run(s))...")
    for repeat in range(1, WARMUPS + 1):
        timed_run(bwa, "warmup_startup", repeat, empty1, empty2)
        timed_run(bwa, "warmup_alignment", repeat, R1, R2)

    rows = []
    print(f"Measuring startup and alignment ({REPEATS} repeat(s), {THREADS} thread(s))...")
    for repeat in range(1, REPEATS + 1):
        startup = timed_run(bwa, "startup", repeat, empty1, empty2)
        rows.append(("startup", repeat, startup))
        alignment = timed_run(bwa, "alignment", repeat, R1, R2)
        rows.append(("alignment", repeat, alignment))
        print(f"  repeat {repeat}: startup {startup:.2f} s, 10k pairs {alignment:.2f} s")

    with (RESULTS_DIR / "timings.tsv").open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("mode", "repeat", "wall_seconds"))
        writer.writerows((mode, repeat, f"{seconds:.6f}") for mode, repeat, seconds in rows)
    return rows


def summarize(rows):
    startup = statistics.median(seconds for mode, _, seconds in rows if mode == "startup")
    batch = statistics.median(seconds for mode, _, seconds in rows if mode == "alignment")
    alignment_only = max(0.0, batch - startup)
    batches = TOTAL_PAIRS // PAIRS
    separate_hours = batches * batch / 3600
    persistent_hours = (startup + batches * alignment_only) / 3600

    summary = [
        ("batch_pairs", PAIRS, "pairs"),
        ("total_pairs", TOTAL_PAIRS, "pairs"),
        ("batches", batches, "invocations"),
        ("threads", THREADS, "threads_per_invocation"),
        ("median_startup", f"{startup:.6f}", "seconds"),
        ("median_10k_batch_total", f"{batch:.6f}", "seconds"),
        ("median_10k_alignment_increment", f"{alignment_only:.6f}", "seconds"),
        ("projected_separate_batches", f"{separate_hours:.6f}", "hours"),
        ("projected_persistent_process", f"{persistent_hours:.6f}", "hours"),
    ]
    with (RESULTS_DIR / "summary.tsv").open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value", "unit"))
        writer.writerows(summary)

    print("\nSummary")
    for metric, value, unit in summary:
        print(f"{metric:35} {value} {unit}")


def main():
    bwa = validate_inputs()
    generate_reads()
    summarize(benchmark(bwa))


if __name__ == "__main__":
    main()
