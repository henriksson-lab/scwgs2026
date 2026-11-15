#!/usr/bin/env python3
"""Download pinned Yanai bacterial references and combine matching FASTA/GFF."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.parse import unquote


def digest(path, algorithm="md5"):
    result = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            result.update(chunk)
    return result.hexdigest()


def publish(temp, output):
    if output.exists():
        if digest(temp, "sha256") != digest(output, "sha256"):
            raise RuntimeError(f"Refusing to replace different reference: {output}")
        temp.unlink()
    else:
        os.link(temp, output)
        temp.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/husky/henriksson/atrandi/bwamem_ref/ling"))
    args = parser.parse_args()
    root = args.root.resolve()
    source = root / "source"
    source.mkdir(parents=True, exist_ok=True)
    here = Path(__file__).resolve().parent
    refs = json.loads((here / "assemblies.json").read_text())
    with (here / "source_manifest.tsv").open() as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    for row in manifest:
        path = source / row["filename"]
        if not path.exists():
            partial = path.with_name(path.name + ".part")
            subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error",
                "--retry", "5", "--continue-at", "-", "--output", str(partial), row["url"]], check=True)
            if partial.stat().st_size != int(row["bytes"]) or digest(partial) != row["md5"]:
                raise RuntimeError(f"NCBI checksum mismatch: {partial}")
            os.link(partial, path)
            partial.unlink()
        if path.stat().st_size != int(row["bytes"]) or digest(path) != row["md5"]:
            raise RuntimeError(f"NCBI checksum mismatch: {path}")
        print("Verified", row["filename"], flush=True)

    sequences = {}
    sequence_species = {}
    with tempfile.NamedTemporaryFile("wb", dir=root, delete=False) as handle:
        fasta_temp = Path(handle.name)
        for ref in refs:
            path = source / (ref["stem"] + "_genomic.fna.gz")
            with gzip.open(path, "rb") as stream:
                current = None
                for line in stream:
                    if line.startswith(b">"):
                        current = line[1:].split()[0].decode()
                        if current in sequences:
                            raise RuntimeError(f"Duplicate FASTA sequence ID {current}")
                        sequences[current] = 0
                        sequence_species[current] = ref["key"]
                    else:
                        if current is None:
                            raise RuntimeError("FASTA sequence before header")
                        sequences[current] += len(line.strip())
                    handle.write(line)
        handle.flush()
        os.fsync(handle.fileno())

    features = {}
    gene_ids = set()
    gff_ids = set()
    gene_names = 0
    with tempfile.NamedTemporaryFile("w", dir=root, delete=False) as handle:
        gff_temp = Path(handle.name)
        handle.write("##gff-version 3\n")
        for ref in refs:
            handle.write(f"# source_assembly {ref['accession']} {ref['species']}\n")
            with gzip.open(source / (ref["stem"] + "_genomic.gff.gz"), "rt") as stream:
                for line in stream:
                    if line.startswith("##gff-version"):
                        continue
                    if line.startswith("##FASTA"):
                        raise RuntimeError("Unexpected embedded FASTA in GFF")
                    if line.startswith("#"):
                        handle.write(line)
                        continue
                    columns = line.rstrip("\n").split("\t")
                    if len(columns) != 9:
                        raise RuntimeError("Malformed GFF record")
                    seqid, _, kind, start, end, *_ = columns
                    if seqid not in sequences or not 1 <= int(start) <= int(end) <= sequences[seqid]:
                        raise RuntimeError(f"GFF sequence/coordinate mismatch: {line}")
                    features[kind] = features.get(kind, 0) + 1
                    if kind == "gene":
                        attributes = dict(part.split("=", 1) for part in columns[8].split(";") if "=" in part)
                        locus = unquote(attributes.get("locus_tag", ""))
                        identifier = unquote(attributes.get("ID", ""))
                        if not locus or locus in gene_ids or not identifier or identifier in gff_ids:
                            raise RuntimeError(f"Missing/duplicate gene identifier: {line}")
                        gene_ids.add(locus)
                        gff_ids.add(identifier)
                        gene_names += bool(attributes.get("Name"))
                    handle.write(line)
        handle.flush()
        os.fsync(handle.fileno())

    publish(fasta_temp, root / "all.fa")
    publish(gff_temp, root / "all.gff")
    shutil.copyfile(here / "assemblies.json", root / "assemblies.json")
    shutil.copyfile(here / "source_manifest.tsv", root / "source_manifest.tsv")
    report = {
        "paper_doi": "10.1038/s41586-023-06974-w",
        "assemblies": refs,
        "contigs": [{"id": name, "length": length, "species": sequence_species[name]}
                    for name, length in sequences.items()],
        "total_bases": sum(sequences.values()),
        "feature_counts": features,
        "unique_gene_locus_tags": len(gene_ids),
        "genes_with_Name": gene_names,
        "bascet_countfeature": {"useFeature": "gene", "attrGeneId": "locus_tag", "attrGeneName": "Name"},
        "all_fa_sha256": digest(root / "all.fa", "sha256"),
        "all_gff_sha256": digest(root / "all.gff", "sha256"),
    }
    (root / "reference.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Ready: {root / 'all.fa'}; {len(sequences)} contigs, {sum(sequences.values())} bases")
    print(f"Annotation: {root / 'all.gff'}; {len(gene_ids)} unique gene locus tags")


if __name__ == "__main__":
    main()
