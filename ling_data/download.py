#!/usr/bin/env python3
"""Download the four Ling/Aarestrup sequencing deposits with Zenodo MD5 checks."""
import argparse
import concurrent.futures
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def md5(path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(row, root):
    folder = root / "source" / row["library"]
    folder.mkdir(parents=True, exist_ok=True)
    output = folder / row["filename"]
    partial = output.with_name(output.name + ".part")
    expected_size = int(row["bytes"])
    if not output.exists():
        print(f"Download {row['library']} A={row['barcode_a']} {row['role']}", flush=True)
        subprocess.run([
            "curl", "--fail", "--location", "--silent", "--show-error",
            "--retry", "8", "--retry-delay", "5", "--connect-timeout", "30",
            "--continue-at", "-", "--output", str(partial), row["url"],
        ], check=True)
        if partial.stat().st_size != expected_size or md5(partial) != row["md5"]:
            raise RuntimeError(f"Zenodo size/MD5 mismatch: {partial}")
        os.link(partial, output)
        partial.unlink()
    if output.stat().st_size != expected_size or md5(output) != row["md5"]:
        raise RuntimeError(f"Zenodo size/MD5 mismatch: {output}")
    receipt = output.with_name(output.name + ".verified.json")
    receipt.write_text(json.dumps({
        "zenodo_record": row["record"], "source_url": row["url"],
        "bytes": expected_size, "md5": row["md5"], "verified_at": time.time(),
        "library": row["library"], "barcode_a": row["barcode_a"], "role": row["role"],
    }, indent=2) + "\n")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/husky/henriksson/atrandi/v9_ling"))
    parser.add_argument("--manifest", type=Path, default=Path(__file__).with_name("manifest.tsv"))
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--library", action="append", help="Limit to a library, repeatable")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".download.lock").open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error(f"Another download is running in {root}")
        rows = list(csv.DictReader(args.manifest.open(), delimiter="\t"))
        if args.library:
            rows = [row for row in rows if row["library"] in args.library]
        if not rows:
            parser.error("No source files selected")
        shutil.copyfile(args.manifest, root / "manifest.tsv")
        shutil.copyfile(Path(__file__).with_name("README.md"), root / "README.md")
        for path in Path(__file__).parent.glob("zenodo_*.json"):
            shutil.copyfile(path, root / path.name)
        print(f"Preparing {len(rows)} source files ({sum(int(r['bytes']) for r in rows)/1e9:.3f} GB)", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(fetch, row, root) for row in rows]
            for count, future in enumerate(concurrent.futures.as_completed(futures), 1):
                print(f"[{count}/{len(futures)}] Verified {future.result().name}", flush=True)
        print(f"All {len(rows)} source FASTQs verified in {root / 'source'}", flush=True)


if __name__ == "__main__":
    main()
