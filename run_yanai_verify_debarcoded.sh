#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
target_dir="${YANAI_VERIFIER_TARGET_DIR:-/tmp/yanai-verifier-target}"
cargo build --locked --offline --release \
  --manifest-path "$repo_dir/yanai_petriseq/verifier/Cargo.toml" --target-dir "$target_dir"
exec "$target_dir/release/verify-yanai-debarcoding" "$@"
