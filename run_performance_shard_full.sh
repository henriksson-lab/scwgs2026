#!/usr/bin/env bash
# Build, benchmark all configured inputs, make figures, and validate the run.
set -euo pipefail

repo_dir="${PERFORMANCE_SHARD_REPO:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
export PERFORMANCE_SHARD_ROOT="${PERFORMANCE_SHARD_ROOT:-/husky/henriksson/atrandi/v9_performance_shard}"
export BASCET_DIR="${BASCET_DIR:-/home/mahogny/github/bascet}"
export ZORN_DIR="${ZORN_DIR:-/home/mahogny/github/zorn}"
export BENCHMARK_OUT="${BENCHMARK_OUT:-$PERFORMANCE_SHARD_ROOT/benchmarks/full_$(date -u +%Y%m%dT%H%M%SZ)}"
export BENCHMARK_ZORN_MODE="${BENCHMARK_ZORN_MODE:-installed}"
export BENCHMARK_THREADS="${BENCHMARK_THREADS:-40}"
export BENCHMARK_MEM="${BENCHMARK_MEM:-50 GB}"
export BENCHMARK_RSS_INTERVAL="${BENCHMARK_RSS_INTERVAL:-0.1}"
export BENCHMARK_NOTE="${BENCHMARK_NOTE:-Full performance run on a quiet system, including fresh repetitions of all earlier test sizes.}"

mkdir -p -- "$(dirname -- "$BENCHMARK_OUT")" "$PERFORMANCE_SHARD_ROOT/benchmarks"
exec 9>"$PERFORMANCE_SHARD_ROOT/benchmarks/.performance-benchmark.lock"
flock -n 9 || { echo 'Another performance benchmark is running' >&2; exit 1; }
if [[ -e "$BENCHMARK_OUT" || -e "$BENCHMARK_OUT.scripts" ]]; then
    echo "Run output already exists: $BENCHMARK_OUT" >&2
    exit 1
fi
trap 'code=$?; if (( code == 0 )); then printf "complete\n"; else printf "failed (exit %s)\n" "$code"; fi > "$BENCHMARK_OUT.status"' EXIT
printf 'building\n' > "$BENCHMARK_OUT.status"
printf '%s\n' "$$" > "$BENCHMARK_OUT.pid"
scripts="$BENCHMARK_OUT.scripts"
mkdir -- "$scripts"
cp -- "$repo_dir/run_performance_shard_benchmark.R" "$scripts/"
cp -- "$repo_dir/supnote_speed_performance_shard.R" "$scripts/"
cp -- "$repo_dir/supnote_speed_performance_shard_rss.R" "$scripts/"
cp -- "${BASH_SOURCE[0]}" "$scripts/driver.sh"
{
    printf 'repository\tcommit\n'
    printf 'bascet\t%s\n' "$(git -C "$BASCET_DIR" rev-parse HEAD)"
    printf 'zorn\t%s\n' "$(git -C "$ZORN_DIR" rev-parse HEAD)"
} > "$scripts/source_commits.tsv"
(
    cd -- "$BASCET_DIR"
    cargo build --release --offline -p bascet-cli --bin bascet -j "${BENCHMARK_BUILD_JOBS:-4}"
)
printf 'verifying inputs and benchmarking\n' > "$BENCHMARK_OUT.status"
Rscript "$scripts/run_performance_shard_benchmark.R"
printf 'plotting\n' > "$BENCHMARK_OUT.status"
Rscript "$scripts/supnote_speed_performance_shard.R" "$BENCHMARK_OUT"
Rscript "$scripts/supnote_speed_performance_shard_rss.R" "$BENCHMARK_OUT"
printf 'validating results\n' > "$BENCHMARK_OUT.status"
Rscript - "$BENCHMARK_OUT" <<'RS'
output <- commandArgs(trailingOnly = TRUE)[1]
meta <- jsonlite::fromJSON(file.path(output, "metadata.json"))
times <- read.delim(file.path(output, "timings.tsv"), check.names = FALSE)
stopifnot(meta$status == "complete", nrow(times) == 2 * sum(meta$plan$repetitions),
          all(times$exit_code == 0), all(is.finite(times$wall_seconds)),
          all(times$wall_seconds >= 0), all(times$peak_rss_bytes > 0))
for(i in seq_len(nrow(meta$plan))) {
  n <- meta$plan$cells[i]
  for(rep in seq_len(meta$plan$repetitions[i])) {
    d <- times[times$cells == n & times[["repeat"]] == rep, ]
    stopifnot(nrow(d) == 2, setequal(d$step, c("debarcode", "shardify")),
              all(d$read_pairs == meta$plan$read_pairs[i]))
    folder <- file.path(output, paste0("cells_", n), paste0("rep_", rep))
    log <- readLines(file.path(folder, "logs", "Zshardify.0.log"), warn = FALSE)
    stopifnot(any(grepl("output_files_created=10", log, fixed = TRUE)))
    if(!isTRUE(meta$keep_outputs)) stopifnot(!dir.exists(file.path(folder, "work")))
  }
}
names <- c("plot_perf_bascet_linear", "plot_perf_bascet_logx",
           "plot_perf_bascet_rss", "plot_perf_bascet_peak_rss")
figures <- file.path(output, "plots", as.vector(outer(names, c(".pdf", ".png"), paste0)))
stopifnot(all(file.exists(figures)), all(file.info(figures)$size > 100))
message <- sprintf("Verified %d successful steps over %d repetitions and %d cell counts; figures complete.",
                   nrow(times), sum(meta$plan$repetitions), nrow(meta$plan))
writeLines(message, file.path(output, "validation.txt"))
cat(message, "\n")
RS
