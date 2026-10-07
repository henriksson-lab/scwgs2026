#!/usr/bin/env bash
set -euo pipefail

ACCESSION="${ACCESSION:-SRR31758484}"
CLEANBAR_DIR="${CLEANBAR_DIR:-/home/mahogny/github/scmetag-revision/supnote_speed/cleanbar/cleanbar}"
CB_BIN="${CB_BIN:-${CLEANBAR_DIR}/CB}"
BARCODES="${BARCODES:-${CLEANBAR_DIR}/barcodes.txt}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="${RUN_DIR:-${SCRIPT_DIR}/runs/cleanbar_${ACCESSION}}"
FIFO_FASTQ="${FIFO_FASTQ:-${RUN_DIR}/${ACCESSION}.fastq}"
SCREEN_LOG="${SCREEN_LOG:-${RUN_DIR}/cleanbar.screen.log}"
TIME_FILE="${TIME_FILE:-${RUN_DIR}/time.txt}"
MANIFEST="${MANIFEST:-${RUN_DIR}/manifest.txt}"
OUTPUT_INVENTORY="${OUTPUT_INVENTORY:-${RUN_DIR}/output_inventory.txt}"
KEEP_OUTPUT="${KEEP_OUTPUT:-0}"

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

if [[ ! -x "${CB_BIN}" ]]; then
    echo "error: CleanBar binary not found or not executable: ${CB_BIN}" >&2
    exit 1
fi

if [[ ! -f "${BARCODES}" ]]; then
    echo "error: CleanBar barcode file not found: ${BARCODES}" >&2
    exit 1
fi

mkdir -p "${RUN_DIR}"
rm -rf "${RUN_DIR}/res_4barcodes" "${RUN_DIR}/res_23barcodes"
rm -f "${FIFO_FASTQ}" \
    "${RUN_DIR}/${ACCESSION}_stats.txt" \
    "${RUN_DIR}/${ACCESSION}_summary.txt" \
    "${RUN_DIR}/${ACCESSION}_1_0_bar.fq" \
    "${RUN_DIR}/${ACCESSION}_3_2_bar.fq" \
    "${SCREEN_LOG}" \
    "${TIME_FILE}" \
    "${OUTPUT_INVENTORY}"

ln -sfn "${CB_BIN}" "${RUN_DIR}/CB"
ln -sfn "${BARCODES}" "${RUN_DIR}/barcodes.txt"
mkdir -p "${RUN_DIR}/res_4barcodes" "${RUN_DIR}/res_23barcodes"
mkfifo "${FIFO_FASTQ}"

cleanup() {
    local status=$?
    if [[ -n "${DECOMPRESS_PID:-}" ]]; then
        wait "${DECOMPRESS_PID}" 2>/dev/null || true
    fi
    rm -f "${FIFO_FASTQ}"
    if [[ "${KEEP_OUTPUT}" != "1" ]]; then
        rm -rf "${RUN_DIR}/res_4barcodes" "${RUN_DIR}/res_23barcodes"
        rm -f "${RUN_DIR}/${ACCESSION}_stats.txt" \
            "${RUN_DIR}/${ACCESSION}_summary.txt" \
            "${RUN_DIR}/${ACCESSION}_1_0_bar.fq" \
            "${RUN_DIR}/${ACCESSION}_3_2_bar.fq"
    fi
    exit "${status}"
}
trap cleanup EXIT INT TERM

cat >"${MANIFEST}" <<MANIFEST_CONTENT
accession=${ACCESSION}
fastq=${FASTQ}
cleanbar_dir=${CLEANBAR_DIR}
cb_bin=${CB_BIN}
barcodes=${BARCODES}
run_dir=${RUN_DIR}
screen_log=${SCREEN_LOG}
time_file=${TIME_FILE}
keep_output=${KEEP_OUTPUT}
note=CleanBar has no thread option; this benchmarks the CB executable as provided.
MANIFEST_CONTENT

echo "Input FASTQ: ${FASTQ}"
echo "Run dir: ${RUN_DIR}"
echo "Screen log: ${SCREEN_LOG}"
echo "Timing log: ${TIME_FILE}"

if [[ "${FASTQ}" == *.gz ]]; then
    gzip -dc "${FASTQ}" >"${FIFO_FASTQ}" &
else
    cat "${FASTQ}" >"${FIFO_FASTQ}" &
fi
DECOMPRESS_PID=$!

(
    cd "${RUN_DIR}"
    set +e
    /usr/bin/time -v -o "${TIME_FILE}" ./CB barcodes.txt "$(basename "${FIFO_FASTQ}")" >"${SCREEN_LOG}" 2>&1
    CB_STATUS=$?
    set -e
    echo "${CB_STATUS}" >cleanbar.exit_status
)

wait "${DECOMPRESS_PID}"
unset DECOMPRESS_PID

CB_STATUS="$(cat "${RUN_DIR}/cleanbar.exit_status")"
{
    echo "cleanbar_exit_status=${CB_STATUS}"
    echo "res_4barcodes_files=$(find "${RUN_DIR}/res_4barcodes" -maxdepth 1 -type f -name '*.fq' | wc -l)"
    echo "res_23barcodes_files=$(find "${RUN_DIR}/res_23barcodes" -maxdepth 1 -type f -name '*.fq' | wc -l)"
    find "${RUN_DIR}" -maxdepth 1 -type f \( -name "${ACCESSION}_*.txt" -o -name "${ACCESSION}_*.fq" \) -print
    du -sh "${RUN_DIR}" 2>/dev/null || true
} >"${OUTPUT_INVENTORY}"

if [[ "${CB_STATUS}" != "0" ]] && ! grep -q "That's all Folks" "${SCREEN_LOG}"; then
    echo "error: CleanBar exited with status ${CB_STATUS}; see ${SCREEN_LOG}" >&2
    exit "${CB_STATUS}"
fi

echo "Done."
echo "Timing log: ${TIME_FILE}"
echo "Screen log: ${SCREEN_LOG}"
echo "Output inventory: ${OUTPUT_INVENTORY}"
if [[ "${KEEP_OUTPUT}" == "1" ]]; then
    echo "Kept CleanBar generated files in ${RUN_DIR}"
else
    echo "Removed bulky CleanBar generated outputs; set KEEP_OUTPUT=1 to retain them."
fi
