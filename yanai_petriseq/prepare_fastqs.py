#!/usr/bin/env python3
"""Verify Yanai ENA FASTQs and stage BGZF mates under GEO library names."""
import argparse
import concurrent.futures
import csv
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time


def checksum(path):
    digest = hashlib.md5()
    with path.open("rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def prepare_mate(row, mate, args):
    run, library = row["run_accession"], row["library_prefix"]
    assert re.fullmatch(r"SRR[0-9]+", run), run
    assert re.fullmatch(r"D[0-9]+_lib[0-9]+", library), library
    urls = row["fastq_ftp"].split(";")
    checksums = row["fastq_md5"].split(";")
    sizes = list(map(int, row["fastq_bytes"].split(";")))
    assert len(urls) == len(checksums) == len(sizes) == 2, run
    i = mate - 1
    assert urls[i].endswith(f"/{run}_{mate}.fastq.gz"), urls[i]
    source = args.raw_dir / f"{run}_{mate}.fastq.gz"
    output = args.output_dir / f"{library}_S1_L001_R{mate}_001.fastq.gz"
    receipt = output.with_name(output.name + ".json")
    identity = dict(run=run, library=library, mate=mate, source_md5=checksums[i])
    if output.exists() and receipt.exists():
        saved = json.loads(receipt.read_text())
        if all(saved.get(k) == v for k, v in identity.items()):
            if output.stat().st_size == saved.get("output_bytes"):
                subprocess.run(["bgzip", "-t", str(output)], check=True)
                return output, True
        raise RuntimeError(f"Output/receipt disagree: {output}")
    if output.exists() or receipt.exists():
        raise RuntimeError(f"Incomplete existing output (refusing overwrite): {output}")
    if not source.exists():
        if not args.download:
            raise FileNotFoundError(f"{source}; pass --download to retrieve it from ENA")
        partial = source.with_name(source.name + ".part")
        print(f"Downloading {run} mate {mate}", flush=True)
        subprocess.run([
            "curl", "--fail", "--location", "--silent", "--show-error",
            "--retry", "8", "--retry-delay", "5", "--connect-timeout", "30",
            "--continue-at", "-", "--output", str(partial), "https://" + urls[i],
        ], check=True)
        if partial.stat().st_size != sizes[i] or checksum(partial) != checksums[i]:
            raise RuntimeError(f"ENA size/MD5 mismatch: {partial}")
        os.rename(partial, source)
    if source.stat().st_size != sizes[i] or checksum(source) != checksums[i]:
        raise RuntimeError(f"ENA size/MD5 mismatch: {source}")

    print(f"BGZF {library} mate {mate}", flush=True)
    fd, name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".part", dir=args.output_dir)
    temp = Path(name)
    try:
        # Inflate the original stream with CRC validation and feed its bytes
        # unchanged to bgzip. Sequences, qualities, read names and order survive.
        with os.fdopen(fd, "wb") as packed:
            proc = subprocess.Popen(
                ["bgzip", "-@", str(args.threads), "-c"], stdin=subprocess.PIPE, stdout=packed
            )
            try:
                with gzip.open(source, "rb") as uncompressed:
                    shutil.copyfileobj(uncompressed, proc.stdin, length=1024 * 1024)
                proc.stdin.close()
                status = proc.wait()
                if status:
                    raise subprocess.CalledProcessError(status, proc.args)
            except BaseException:
                if proc.poll() is None:
                    proc.terminate()
                proc.wait()
                raise
            packed.flush()
            os.fsync(packed.fileno())
        subprocess.run(["bgzip", "-t", str(temp)], check=True)
        # Publish without replacing an output produced by another invocation.
        os.link(temp, output)
        temp.unlink()
        identity.update(output_bytes=output.stat().st_size, completed_at=time.time())
        with tempfile.NamedTemporaryFile("w", dir=args.output_dir, delete=False) as handle:
            json.dump(identity, handle, indent=2)
            handle.write("\n")
            receipt_temp = Path(handle.name)
        os.link(receipt_temp, receipt)
        receipt_temp.unlink()
        print(f"Ready {output.name}", flush=True)
        return output, False
    finally:
        temp.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path(__file__).with_name("fastq_manifest.tsv"))
    parser.add_argument("--root", type=Path, required=True, help="Dataset root; creates raw_gzip and raw_bgzf")
    parser.add_argument("--raw-dir", type=Path, help="Existing native SRR FASTQ directory")
    parser.add_argument("--output-dir", type=Path, help="Override directory for Illumina-style BGZF names")
    parser.add_argument("--download", action="store_true", help="Download missing mates, checking ENA MD5")
    parser.add_argument("--workers", type=int, default=4, help="Concurrent FASTQ downloads/recompressions")
    parser.add_argument("--threads", type=int, default=2, help="bgzip threads per worker")
    parser.add_argument("--run", action="append", dest="runs", help="Restrict to a run; repeatable")
    args = parser.parse_args()
    if args.workers < 1 or args.threads < 1:
        parser.error("--workers and --threads must be positive")
    args.root = args.root.resolve()
    args.raw_dir = (args.raw_dir or args.root / "raw_gzip").resolve()
    args.output_dir = (args.output_dir or args.root / "raw_bgzf").resolve()
    for directory in [args.root, args.raw_dir, args.output_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    # Protect resumable .part files against concurrent invocations. The PID
    # check also covers a downloader started before lock support was added.
    pidfile = args.root / "download.pid"
    if pidfile.exists():
        try:
            pid = int(pidfile.read_text().strip())
            running = (Path("/proc") / str(pid) / "cmdline").read_bytes().split(b"\0")
            if pid != os.getpid() and str(Path(__file__).resolve()).encode() in running \
                    and str(args.root).encode() in running:
                parser.error(f"Another Yanai downloader is running (PID {pid})")
        except (OSError, ValueError):
            pass
    lock = (args.root / ".download.lock").open("w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        parser.error(f"Another Yanai download is running in {args.root}")
    pidfile.write_text(str(os.getpid()) + "\n")
    rows = list(csv.DictReader(args.manifest.open(), delimiter="\t"))
    if args.runs:
        requested = set(args.runs)
        rows = [r for r in rows if r["run_accession"] in requested]
        if len(rows) != len(requested):
            parser.error("Some requested --run accessions are absent from the manifest")
    if not rows:
        parser.error("No runs selected")
    if len({r["library_prefix"] for r in rows}) != len(rows):
        parser.error("Manifest must have a unique library prefix for each run")
    (args.root / "rawdata.txt").write_text(str(args.output_dir) + "\n")
    if args.manifest.resolve() != (args.root / "fastq_manifest.tsv").resolve():
        shutil.copyfile(args.manifest, args.root / "fastq_manifest.tsv")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(prepare_mate, row, mate, args) for row in rows for mate in [1, 2]]
        complete = 0
        for future in concurrent.futures.as_completed(futures):
            path, reused = future.result()
            complete += 1
            print(f"[{complete}/{len(futures)}] {'Verified existing' if reused else 'Prepared'} {path.name}", flush=True)
    print(f"All {len(rows)} libraries ready; Zorn rawdata.txt -> {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
