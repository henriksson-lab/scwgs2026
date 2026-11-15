#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
build_dir="${LING_CARGO_TARGET_DIR:-/tmp/ling2tirp-target}"
cargo build --locked --release --manifest-path "$repo_dir/ling2tirp/Cargo.toml" \
  --target-dir "$build_dir"
exec "$build_dir/release/ling2tirp" --manifest "$repo_dir/ling_data/manifest.tsv" "$@"
