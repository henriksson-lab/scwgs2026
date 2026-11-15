#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
ref_dir="/husky/henriksson/atrandi/bwamem_ref/ling"
if [[ $# -gt 0 ]]; then ref_dir="$1"; shift; fi
if [[ $# -gt 0 ]]; then
    printf 'Usage: %s [reference_directory]\n' "$0" >&2
    exit 2
fi
bascet_bin="/home/mahogny/github/bascet/target/release/bascet"
if [[ -n "$(printenv BASCET_BIN || true)" ]]; then bascet_bin="$(printenv BASCET_BIN)"; fi
test -s "$ref_dir/all.fa"
index_count=0
for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
    if [[ -s "$ref_dir/all.fa$suffix" ]]; then index_count=$((index_count + 1)); fi
done
if [[ "$index_count" -eq 5 ]]; then
    printf 'Existing BWA-MEM2 index found: %s/all.fa\n' "$ref_dir"
elif [[ "$index_count" -ne 0 ]]; then
    printf 'Incomplete existing index; refusing to overwrite it: %s\n' "$ref_dir" >&2
    exit 1
else
    "$bascet_bin" exttool bwa-mem2 index "$ref_dir/all.fa" 2>&1 | tee "$ref_dir/index.log"
fi
for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do test -s "$ref_dir/all.fa$suffix"; done
Rscript "$repo_dir/yanai_reference/record_index.R" "$ref_dir" "$bascet_bin"
printf 'Bascet-compatible BWA-MEM2 index ready: %s/all.fa\n' "$ref_dir"
