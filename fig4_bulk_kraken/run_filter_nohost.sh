#!/usr/bin/env bash
set -euo pipefail

base="/husky/henriksson/atrandi/rawdata/251128_jyoti_mock_bulk/bulk/after_filter"

cargo run --release -- filterfastq \
  -@ 20 \
  --keep unmapped \
  --in "${base}/bulk_unsorted.sam" \
  --out-r1 "${base}/filtered_bulk.R1.fq.gz" \
  --out-r2 "${base}/filtered_bulk.R2.fq.gz"
