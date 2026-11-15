#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export GORG_ROOT="${GORG_ROOT:-/husky/henriksson/atrandi/v9_wgs_gorg}"
export GORG_RUNNER=local
export GORG_THREADS="${GORG_THREADS:-12}"
export GORG_MEM="${GORG_MEM:-24g}"
exec 9>"${GORG_ROOT}/.all-workflow.lock"
flock -n 9 || { echo 'Another GORG workflow is running' >&2; exit 1; }
cd "$repo_dir"
Rscript prep_gorg.R
exec 8>"/husky/henriksson/atrandi/.umap-analysis.lock"
flock 8
Rscript post_gorg.R
Rscript verify_umaps.R gorg
