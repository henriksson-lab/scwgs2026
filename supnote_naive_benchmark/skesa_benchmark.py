#!/usr/bin/env python3
"""Benchmark SKESA startup and assembly of the same 10k read pairs."""

import csv
import shutil
import statistics
import subprocess
import time
from pathlib import Path


# Keep these values in sync with bwamem2_benchmark.py.
SKESA = "skesa"
PAIRS = 10_000
TOTAL_PAIRS = 800_000_000
THREADS = 10
REPEATS = 5
WARMUPS = 1

R1 = Path(f"data/reads_{PAIRS}_R1.fastq.gz")
R2 = Path(f"data/reads_{PAIRS}_R2.fastq.gz")
RESULTS_DIR = Path("results")


def timed_run(command, log):
    start = time.perf_counter()
    with log.open("w") as stderr:
        completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=stderr)
    elapsed = time.perf_counter() - start
    if completed.returncode:
        raise RuntimeError(f"SKESA failed; see {log}")
    return elapsed


def main():
    skesa = shutil.which(SKESA)
    if skesa is None:
        raise RuntimeError(f"cannot find {SKESA} on PATH")
    if not R1.is_file() or not R2.is_file():
        from bwamem2_benchmark import generate_reads

        generate_reads()
    if TOTAL_PAIRS % PAIRS:
        raise RuntimeError("TOTAL_PAIRS must be divisible by PAIRS")

    logs = RESULTS_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    startup_command = [skesa, "--version"]
    assembly_command = [
        skesa,
        "--reads",
        f"{R1},{R2}",
        "--cores",
        str(THREADS),
        "--contigs_out",
        "/dev/null",
    ]

    print(f"Warming up ({WARMUPS} run(s))...")
    for repeat in range(1, WARMUPS + 1):
        timed_run(startup_command, logs / f"skesa_warmup_startup_{repeat}.log")
        timed_run(assembly_command, logs / f"skesa_warmup_assembly_{repeat}.log")

    rows = []
    print(f"Measuring SKESA ({REPEATS} repeat(s), {THREADS} thread(s))...")
    for repeat in range(1, REPEATS + 1):
        startup = timed_run(startup_command, logs / f"skesa_startup_{repeat}.log")
        assembly = timed_run(assembly_command, logs / f"skesa_assembly_{repeat}.log")
        rows.extend((("startup", repeat, startup), ("assembly", repeat, assembly)))
        print(f"  repeat {repeat}: startup {startup:.4f} s, 10k pairs {assembly:.2f} s")

    timings_path = RESULTS_DIR / "skesa_timings.tsv"
    with timings_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("mode", "repeat", "wall_seconds"))
        writer.writerows((mode, repeat, f"{seconds:.6f}") for mode, repeat, seconds in rows)

    startup = statistics.median(seconds for mode, _, seconds in rows if mode == "startup")
    batch = statistics.median(seconds for mode, _, seconds in rows if mode == "assembly")
    assembly_only = max(0.0, batch - startup)
    batches = TOTAL_PAIRS // PAIRS
    separate_hours = batches * batch / 3600
    persistent_hours = (startup + batches * assembly_only) / 3600
    summary = [
        ("batch_pairs", PAIRS, "pairs"),
        ("total_pairs", TOTAL_PAIRS, "pairs"),
        ("batches", batches, "invocations"),
        ("threads", THREADS, "threads_per_invocation"),
        ("median_startup", f"{startup:.6f}", "seconds"),
        ("median_10k_batch_total", f"{batch:.6f}", "seconds"),
        ("median_10k_assembly_increment", f"{assembly_only:.6f}", "seconds"),
        ("projected_separate_batches", f"{separate_hours:.6f}", "hours"),
        ("projected_persistent_process", f"{persistent_hours:.6f}", "hours"),
    ]
    summary_path = RESULTS_DIR / "skesa_summary.tsv"
    with summary_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value", "unit"))
        writer.writerows(summary)

    print("\nSummary")
    for metric, value, unit in summary:
        print(f"{metric:35} {value} {unit}")
    print(f"\nWrote {timings_path} and {summary_path}")


if __name__ == "__main__":
    main()
