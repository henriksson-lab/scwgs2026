#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
download_root="${YANAI_ROOT:-/husky/henriksson/atrandi/v9_rnaseq_yanai}"
arguments=("$@")
for ((i=0; i<${#arguments[@]}; i++)); do
    case "${arguments[i]}" in
        --root) download_root="${arguments[i+1]:?Missing --root value}" ;;
        --root=*) download_root="${arguments[i]#--root=}" ;;
    esac
done
mkdir -p "$download_root"
exec flock -n "$download_root/.download.lock" Rscript "${repo_dir}/yanai_petriseq/prepare_fastqs.R" \
    --root "$download_root" \
    --workers "${YANAI_DOWNLOAD_WORKERS:-8}" \
    --threads "${YANAI_BGZIP_THREADS:-2}" \
    --download "$@"
