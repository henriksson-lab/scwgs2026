#!/usr/bin/env bash
set -euo pipefail

ACCESSION="${ACCESSION:-SRR31758484}"
OUTDIR="${OUTDIR:-data/fastq/${ACCESSION}}"
THREADS="${THREADS:-8}"
TMPDIR="${TMPDIR:-${OUTDIR}/tmp}"

need_cmd() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "error: required command '$1' was not found in PATH" >&2
        exit 1
    fi
}

need_cmd prefetch
need_cmd fasterq-dump

if command -v pigz >/dev/null 2>&1; then
    GZIP_CMD=(pigz -p "${THREADS}")
else
    need_cmd gzip
    GZIP_CMD=(gzip)
fi

mkdir -p "${OUTDIR}" "${TMPDIR}"

if find "${OUTDIR}" -maxdepth 1 -type f -name "${ACCESSION}*.fastq.gz" | grep -q .; then
    echo "FASTQ files already exist in ${OUTDIR}; remove them to download again." >&2
    find "${OUTDIR}" -maxdepth 1 -type f -name "${ACCESSION}*.fastq.gz" -print
    exit 0
fi

echo "Downloading ${ACCESSION} with prefetch..."
prefetch "${ACCESSION}" --output-directory "${OUTDIR}"

echo "Converting ${ACCESSION} to FASTQ..."
fasterq-dump \
    "${OUTDIR}/${ACCESSION}" \
    --outdir "${OUTDIR}" \
    --temp "${TMPDIR}" \
    --split-files \
    --threads "${THREADS}" \
    --progress

echo "Compressing FASTQ files..."
find "${OUTDIR}" -maxdepth 1 -type f -name "${ACCESSION}*.fastq" -print0 |
    xargs -0 -r -n 1 "${GZIP_CMD[@]}"

echo "Created:"
find "${OUTDIR}" -maxdepth 1 -type f -name "${ACCESSION}*.fastq.gz" -print
