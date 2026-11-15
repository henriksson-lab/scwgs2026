#!/usr/bin/env bash
# Run the existing complete workflow with R temporary files on the dataset disk.
# This changes temporary storage only; all prep/post/audit stages are unchanged.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
dataset="${1:-}"
case "$dataset" in
    gorg) dataset_root="${GORG_ROOT:-/husky/henriksson/atrandi/v9_wgs_gorg}" ;;
    ling) dataset_root="${LING_ROOT:-/husky/henriksson/atrandi/v9_ling}" ;;
    yanai) dataset_root="${YANAI_ROOT:-/husky/henriksson/atrandi/v9_rnaseq_yanai}" ;;
    *) printf 'Usage: %s {gorg|ling|yanai}\n' "$0" >&2; exit 2 ;;
esac
if [[ $# -ne 1 ]]; then
    printf 'Usage: %s {gorg|ling|yanai}\n' "$0" >&2
    exit 2
fi
export TMPDIR="${DATASET_TMPDIR:-${dataset_root}/temp/R_workflow}"
mkdir -p "$TMPDIR"
exec bash "$repo_dir/run_${dataset}_all.sh"
