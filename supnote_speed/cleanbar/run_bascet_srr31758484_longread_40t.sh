#!/usr/bin/env bash
set -euo pipefail

ACCESSION="${ACCESSION:-SRR31758484}"
THREADS="${THREADS:-40}"
MEMORY="${MEMORY:-30GB}"
BASCET_DIR="${BASCET_DIR:-/home/mahogny/github/bascet}"
BASCET_BIN="${BASCET_BIN:-${BASCET_DIR}/target/release/bascet}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="${RUN_DIR:-${SCRIPT_DIR}/runs/bascet_${ACCESSION}_atrandi_wgslr_${THREADS}t}"
TEMP_DIR="${TEMP_DIR:-${RUN_DIR}/tmp}"
OUT_TIRP="${OUT_TIRP:-${RUN_DIR}/${ACCESSION}.atrandi_wgslr.${THREADS}t.tirp.gz}"
LOG_FILE="${LOG_FILE:-${RUN_DIR}/bascet.log}"
TIME_FILE="${TIME_FILE:-${RUN_DIR}/time.txt}"
HIST_FILE="${HIST_FILE:-${RUN_DIR}/${ACCESSION}.hist.tsv}"

if [[ -z "${FASTQ:-}" ]]; then
    for candidate in \
        "${SCRIPT_DIR}/data/fastq/${ACCESSION}/${ACCESSION}.fastq.gz" \
        "${SCRIPT_DIR}/data/fastq/${ACCESSION}/${ACCESSION}_1.fastq.gz" \
        "/husky/henriksson/atrandi/rawdata/cleanbar_longread/${ACCESSION}.fastq.gz"
    do
        if [[ -f "${candidate}" ]]; then
            FASTQ="${candidate}"
            break
        fi
    done
fi

if [[ -z "${FASTQ:-}" || ! -f "${FASTQ}" ]]; then
    echo "error: FASTQ not found. Set FASTQ=/path/to/${ACCESSION}.fastq.gz or run ./download_srr31758484_fastq.sh first." >&2
    exit 1
fi

if [[ ! -x "${BASCET_BIN}" ]]; then
    echo "Bascet release binary not found at ${BASCET_BIN}; building it in ${BASCET_DIR}..." >&2
    cargo build --release --manifest-path "${BASCET_DIR}/Cargo.toml"
fi

mkdir -p "${RUN_DIR}" "${TEMP_DIR}"

cat >"${RUN_DIR}/manifest.txt" <<MANIFEST
accession=${ACCESSION}
fastq=${FASTQ}
bascet_dir=${BASCET_DIR}
bascet_bin=${BASCET_BIN}
threads=${THREADS}
memory=${MEMORY}
temp_dir=${TEMP_DIR}
out_tirp=${OUT_TIRP}
hist_file=${HIST_FILE}
log_file=${LOG_FILE}
time_file=${TIME_FILE}
chemistry=atrandi-wgslr
MANIFEST

echo "Input FASTQ: ${FASTQ}"
echo "Output shard: ${OUT_TIRP}"
echo "Timing log: ${TIME_FILE}"
echo "Bascet log: ${LOG_FILE}"

/usr/bin/time -v -o "${TIME_FILE}" \
    "${BASCET_BIN}" \
    --log-mode="${LOG_FILE}" \
    --log-level=info \
    debarcode \
    --temp="${TEMP_DIR}" \
    --threads="${THREADS}" \
    --memory="${MEMORY}" \
    --r1="${FASTQ}" \
    --out="${OUT_TIRP}" \
    --hist="${HIST_FILE}" \
    atrandi-wgslr

echo "Done."
echo "Output shard: ${OUT_TIRP}"
echo "Histogram: ${HIST_FILE}"
echo "Timing log: ${TIME_FILE}"
