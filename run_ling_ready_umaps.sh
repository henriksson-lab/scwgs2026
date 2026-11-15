#!/usr/bin/env bash
# Run the two ready reductions while the full workflow counts informative k-mers.
# The full workflow later reuses these results and adds the informative UMAP.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export LING_ROOT="${LING_ROOT:-/husky/henriksson/atrandi/v9_ling}"
export LING_UMAP_METHODS=kraken,countsketch
exec 8>"/husky/henriksson/atrandi/.umap-analysis.lock"
flock 8
cd "$repo_dir"
Rscript post_ling.R
