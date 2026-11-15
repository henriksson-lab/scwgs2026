#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
download_root="${LING_ROOT:-/husky/henriksson/atrandi/v9_ling}"
arguments=("$@")
for ((i=0; i<${#arguments[@]}; i++)); do
    case "${arguments[i]}" in
        --root) download_root="${arguments[i+1]:?Missing --root value}" ;;
        --root=*) download_root="${arguments[i]#--root=}" ;;
    esac
done
mkdir -p "$download_root"
exec flock -n "$download_root/.download.lock" Rscript "$repo_dir/ling_data/download.R" \
    --root "$download_root" "$@"
