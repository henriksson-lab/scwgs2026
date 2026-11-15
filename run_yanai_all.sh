#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export YANAI_ROOT="${YANAI_ROOT:-/husky/henriksson/atrandi/v9_rnaseq_yanai}"
export YANAI_RUNNER=local
export YANAI_THREADS="${YANAI_THREADS:-10}"
export YANAI_MEM="${YANAI_MEM:-24g}"
export YANAI_ARRAY_WORKERS="${YANAI_ARRAY_WORKERS:-2}"
exec 9>"${YANAI_ROOT}/.all-workflow.lock"
flock -n 9 || { echo 'Another Yanai workflow is running' >&2; exit 1; }
cd "$repo_dir"
Rscript prep_yanai.R
exec 8>"/husky/henriksson/atrandi/.umap-analysis.lock"
flock 8
Rscript post_yanai.R
Rscript verify_umaps.R yanai
