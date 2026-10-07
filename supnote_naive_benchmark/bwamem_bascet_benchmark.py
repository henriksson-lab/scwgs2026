#!/usr/bin/env python3
"""Benchmark the BWA-MEM2 CLI embedded in the release Bascet binary."""

import csv
import statistics
import subprocess
import time
from pathlib import Path

from bwamem2_benchmark import REFERENCE, generate_reads


# This is intentionally the release binary, not a bwa-mem2 executable on PATH.
BASCET = Path.home() / "github/bascet/target/release/bascet"
PAIRS = 10_000
TOTAL_PAIRS = 800_000_000
THREADS = 10
REPEATS = 5
WARMUPS = 1

R1 = Path(f"data/reads_{PAIRS}_R1.fastq.gz")
R2 = Path(f"data/reads_{PAIRS}_R2.fastq.gz")
RESULTS_DIR = Path("results")
LOGS_DIR = RESULTS_DIR / "logs"


def validate_inputs():
    if not BASCET.is_file():
        raise RuntimeError(f"missing Bascet release binary: {BASCET}")
    if not R1.is_file() or not R2.is_file():
        generate_reads()
    for path in (R1, R2):
        if not path.is_file():
            raise RuntimeError(f"missing input: {path}")
    for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"):
        index = Path(f"{REFERENCE}{suffix}")
        if not index.is_file():
            raise RuntimeError(f"missing bwa-mem2 index: {index}")
    if TOTAL_PAIRS % PAIRS:
        raise RuntimeError("TOTAL_PAIRS must be divisible by PAIRS")


def timed_run(label, repeat, read1, read2):
    command = [
        str(BASCET),
        "exttool",
        "bwa-mem2",
        "mem",
        "-t",
        str(THREADS),
        str(REFERENCE),
        str(read1),
        str(read2),
    ]
    log = LOGS_DIR / f"bascet_bwamem_{label}_{repeat}.log"
    start = time.perf_counter()
    with log.open("w") as stderr:
        completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=stderr)
    elapsed = time.perf_counter() - start
    if completed.returncode:
        raise RuntimeError(f"{label} failed; see {log}")
    return elapsed


def benchmark():
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    empty1 = RESULTS_DIR / "empty_R1.fastq"
    empty2 = RESULTS_DIR / "empty_R2.fastq"
    empty1.write_text("")
    empty2.write_text("")

    print(f"Warming up embedded BWA-MEM2 ({WARMUPS} run(s))...")
    for repeat in range(1, WARMUPS + 1):
        timed_run("warmup_startup", repeat, empty1, empty2)
        timed_run("warmup_alignment", repeat, R1, R2)

    rows = []
    print(f"Measuring embedded BWA-MEM2 ({REPEATS} repeat(s), {THREADS} thread(s))...")
    for repeat in range(1, REPEATS + 1):
        startup = timed_run("startup", repeat, empty1, empty2)
        rows.append(("startup", repeat, startup))
        alignment = timed_run("alignment", repeat, R1, R2)
        rows.append(("alignment", repeat, alignment))
        print(f"  repeat {repeat}: startup {startup:.2f} s, 10k pairs {alignment:.2f} s")

    timings_path = RESULTS_DIR / "bascet_bwamem_timings.tsv"
    with timings_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("mode", "repeat", "wall_seconds"))
        writer.writerows((mode, repeat, f"{seconds:.6f}") for mode, repeat, seconds in rows)
    return rows, timings_path


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
    summary_path = RESULTS_DIR / "bascet_bwamem_summary.tsv"
    with summary_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value", "unit"))
        writer.writerows(summary)

    print("\nSummary")
    for metric, value, unit in summary:
        print(f"{metric:35} {value} {unit}")
    return summary_path


def main():
    validate_inputs()
    rows, timings_path = benchmark()
    summary_path = summarize(rows)
    print(f"\nWrote {timings_path} and {summary_path}")


if __name__ == "__main__":
    main()
