#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
exec Rscript "$repo_dir/yanai_reference/prepare.R" "$@"
