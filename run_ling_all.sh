#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export LING_ROOT="${LING_ROOT:-/husky/henriksson/atrandi/v9_ling}"
export LING_RUNNER=local
export LING_THREADS="${LING_THREADS:-20}"
export LING_MEM="${LING_MEM:-32g}"
mkdir -p "$LING_ROOT"
exec 9>"${LING_ROOT}/.all-workflow.lock"
flock -n 9 || { echo 'Another full Ling workflow is running' >&2; exit 1; }
cd "$repo_dir"
bash run_ling_import.sh --root "$LING_ROOT" \
  --workers "${LING_IMPORT_WORKERS:-8}" --sort-memory "${LING_SORT_MEMORY:-4G}"
Rscript prep_ling.R
exec 8>"/husky/henriksson/atrandi/.umap-analysis.lock"
flock 8
Rscript post_ling.R
Rscript verify_umaps.R ling
