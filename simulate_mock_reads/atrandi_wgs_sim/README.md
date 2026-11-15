# Atrandi WGS sorting benchmark simulator

Build and run from the repository root:

```sh
cargo build --release --locked --manifest-path simulate_mock_reads/atrandi_wgs_sim/Cargo.toml
simulate_mock_reads/atrandi_wgs_sim/target/release/atrandi-wgs-sim \
  --cells 10000 --reads 1000000 --seed 42 --out-prefix /tmp/atrandi_sim
```

This writes `/tmp/atrandi_sim_R1.fastq.gz` and `/tmp/atrandi_sim_R2.fastq.gz`.
`--reads` counts **read pairs**: the example produces one million records in
each FASTQ. Both mates are exactly 150 bases, including the R2 barcode prefix.
Output is BGZF-compressed (using libdeflate), ready for Bascet's debarcoder.
Compression runs in parallel with eight workers per mate by default;
`--compression-threads` sets this worker count. Worker count does not change
the seeded FASTQ contents or their order.
Use `--uncompressed` to write plain `_R1.fastq` and `_R2.fastq` instead.
Existing output files are not overwritten.
Parent output directories are created as needed.

Distinct full barcode combinations are chosen uniformly without replacement.
Each read pair then independently samples one of these cells uniformly, giving
random cell order and approximately equal reads per cell. At low reads/cell,
some selected cells may have no reads; there is no guarantee of observing all
cells. Matching headers contain the read pair index and selected cell index
(`@sim:123:cell42/1` and `/2`); the cell index is only a truth label, not a
Bascet barcode name. A fixed seed reproduces the same output with the same
binary and arguments.

R1 contains random A/C/G/T bases. R2 has this structure:

```text
round4[8] AGGA round3[8] ACTC round2[8] AAGG round1[8] T randomDNA[105]
```

All qualities are `F` (Phred+33 Q37). There are no sequencing errors, UMIs, or
adapter read-through. The genomic portions of the mates are independent random
DNA, representing nonoverlapping ends for a sorting benchmark; they do not
come from a reference genome or model biological insert lengths. Generation
streams through buffered writers, using memory proportional to the number of
cells, not the number of reads.

The embedded `atrandi_barcodes.tsv` was copied from
`/home/mahogny/github/bascet/crates/bascet-cli/src/barcode/atrandi_barcodes.tsv`
at Bascet commit `c378038b1a956c3e7312f680329340c9362fd693`. Barcode order and
positions follow `atrandi_wgs_barcode_illumina.rs`; linker sequences follow
the accompanying `atrandi_wgs_barcode_longread.rs` chemistry comments.
There are 24 barcodes per round, so `--cells` must be in `1..=331776`.
The whitelist is embedded at compile time; running the simulator does not
require the Bascet checkout.

Validation: simulator tests check all 331,776 distinct barcode combinations
and FASTQ record structure. An end-to-end run through the local Bascet
debarcoder recovered the expected 100 cell barcode names from 100,000 input
pairs and retained 99,999 pairs. Independent random genomic sequences can
occasionally trigger Bascet's overlap/trimming filters.

Example Bascet invocation (current local CLI):

```sh
/home/mahogny/github/bascet/target/release/bascet debarcode \
  -1 /tmp/atrandi_sim_R1.fastq.gz -2 /tmp/atrandi_sim_R2.fastq.gz \
  -o /tmp/atrandi_sim.tirp.bbgz --library sim -@ 6 -m '1 GiB' atrandi-wgs
```

Run the simulator's tests with:

```sh
cargo test --locked --manifest-path simulate_mock_reads/atrandi_wgs_sim/Cargo.toml
```

## Performance benchmark datasets

From the repository root, run:

```sh
./run_performance_shard_simulations.sh --dry-run
./run_performance_shard_simulations.sh
```

The script builds the release binary with `Cargo.lock` and generates datasets
under `/husky/henriksson/atrandi/v9_performance_shard/cells_<count>/` for
1, 10, 100, 1,000, 5,000, 10,000, and 50,000 cells. Each cell contributes an expected
11,765 pairs: approximately 1x coverage of a 3 Mb genome using 255 genomic
bases per pair. The datasets contain 11,765; 117,650; 1,176,500; 11,765,000;
58,825,000; 117,650,000; and 588,250,000 read pairs, respectively.

The final full-lane input uses 276,243 cells and exactly 3,250,000,000 read
pairs, representing one Illumina 25B lane at approximately the same coverage
(11,765.004 pairs per cell on average). This dataset uses the exact lane total
rather than rounding the per-cell count before multiplication.

Each dataset has `reads_R1.fastq.gz`, `reads_R2.fastq.gz` (BGZF),
`parameters.tsv`, `simulator.log`, and `checksums.sha256`. Seed 42 is used for
each dataset. Parameters record source and binary hashes and the Rust compiler
version. Completed datasets are published by renaming a temporary directory;
reruns verify parameters and FASTQ checksums before skipping existing datasets.
Mismatched or damaged datasets cause an error and are preserved.

`SIM_COMPRESSION_THREADS` controls workers per mate (default 8), and is recorded
in each dataset's parameters. `PERFORMANCE_SHARD_ROOT`, `SIM_SEED`, and the space-separated `SIM_CELL_COUNTS`
can be overridden, for example for a small local check:

```sh
PERFORMANCE_SHARD_ROOT=/tmp/performance_shard_check SIM_CELL_COUNTS='1 10' \
  ./run_performance_shard_simulations.sh
```

## Timing debarcoding and sharding

For a complete run with automatic plotting and result checks, use:

```sh
./run_performance_shard_full.sh
```

This builds the release binary, archives the exact benchmark/plot scripts and
source commit IDs beside the run, executes all configured sizes, makes both
runtime and RSS figures, and checks every repetition and shard count. A shared
lock prevents two full benchmarks from running simultaneously. Adjacent
`.pid` and `.status` files record the driver PID and its current stage; final
results include `validation.txt`. Intermediate read data are removed according
to the benchmark's normal cleanup policy. Standard benchmark environment
overrides apply, including `BENCHMARK_OUT` for a new output directory.

To run only the benchmark, without the build/plot/validation driver:

```sh
Rscript run_performance_shard_benchmark.R --dry-run
Rscript run_performance_shard_benchmark.R
```

The benchmark uses 40 threads and a 50 GB memory budget for each Bascet step.
For each dataset it runs debarcoding/sorting followed by sharding into exactly
10 files. It uses three repetitions for cell counts up to and including 5,000,
and two for 10,000 and 50,000. The full-lane dataset (276,243 cells,
3,250,000,000 pairs) uses one repetition. Each repetition starts with fresh outputs. Steps and
repetitions run sequentially. A single Bascet binary is used throughout, with
debarcoding compression level 6 and library name `sim`.

Results default to a new timestamped directory under
`/husky/henriksson/atrandi/v9_performance_shard/benchmarks/`. Set `BENCHMARK_OUT`
to choose a new results directory. The R script calls Zorn's `BascetDebarcode`
and `BascetShardify` through an asynchronous `LocalRunner`, supplying the full
50 GB explicitly to both wrappers. Steps still run sequentially. GNU time
records each Bascet subprocess's wall time and peak RSS, excluding the R
launcher/worker. `timings.tsv` records these as `wall_seconds` and
`peak_rss_bytes`; `wrapper_wall_seconds` records the entire async Zorn call
and wait separately. Process startup and I/O are included in Bascet wall time.
Errors are recorded with `exit_code=1` (success is 0) and stop the run.
`summary.tsv` contains successful repetition counts, medians, minima, and maxima.
`metadata.json` records parameters, input hashes, binary hash/version, host,
CPU information, cache policy, and completion status. Zorn source hashes,
`sessionInfo.txt`, wrapper settings (`zorn_parameters.rds`), sharding plans,
and Bascet logs are also saved. All recovered cells are included in sharding,
including the one-cell dataset; no knee/count cutoff is applied.

While each job is running, the R controller samples only its Bascet descendant
PID (not the R worker or individual threads) and writes `rss.tsv`, plus separate
`debarcode.rss.tsv` and `shardify.rss.tsv` traces in each repetition directory.
The requested interval defaults to 0.1 seconds and can be set using
`BENCHMARK_RSS_INTERVAL`. The actual sampling timestamps are recorded.
Very short steps may fall between samples; GNU time's independent peak RSS
measurement still records their peak. Sampling uses ps/procfs RSS, which is a
lightweight approximate resident-memory counter. This profiling requires Linux
and `/usr/bin/time`.

Input checksums are verified before timing; caches are not cleared and no
additional warmup is performed. Sharding follows debarcoding immediately, so
its input may be cached. Generation, verification, and cleanup are excluded
from the reported times, as is preparation of the all-cell sharding plan.
Large intermediate outputs are removed after each successful repetition;
set `BENCHMARK_KEEP_OUTPUTS=1` to retain them. Failed repetition
outputs are kept for diagnosis. Existing result directories are never reused.

For a small smoke run (still three repetitions):

```sh
PERFORMANCE_SHARD_ROOT=/tmp/performance_shard_check BENCHMARK_CELLS=1 \
  BENCHMARK_MEM='2 GB' BENCHMARK_OUT=/tmp/performance_shard_timing_check \
  Rscript run_performance_shard_benchmark.R
```

Other overrides are `BENCHMARK_THREADS` (default 40), `BENCHMARK_MEM` (default
`50 GB`), space-separated `BENCHMARK_CELLS`, `ZORN_DIR` (default
`~/github/zorn`), and `BASCET_DIR` (default `~/github/bascet`, release build).
The default `BENCHMARK_ZORN_MODE=source` loads the Zorn checkout, as in the
repository's other workflows. Set `BENCHMARK_ZORN_MODE=installed` to use
`library(Zorn)` instead. `BENCHMARK_NOTE` records a free-text run annotation.

Plot a completed run (or successful timings from an interrupted run) with:

```sh
Rscript supnote_speed_performance_shard.R /path/to/benchmark_run
Rscript supnote_speed_performance_shard_rss.R /path/to/benchmark_run
```

The R script accepts either the run directory or its `timings.tsv`, plus an
optional output directory. It writes PDF and PNG plots in `<run>/plots/`,
with separate panels for debarcoding/sorting and 10-way sharding. Points show
individual repetitions, diamonds/lines show medians, and ranges show minima
and maxima. Linear-linear, log-cell-axis, and log-log versions are saved,
along with `performance_shard_summary.tsv`. Runtime plots use hours and
independent y scales; the underlying logged timings and summary stay in seconds.
The RSS script produces memory-over-time panels for each cell count and step,
with separate lines for repetitions, and a companion peak-RSS-versus-cell-count
plot. All memory plots are saved as PDF and PNG beside the runtime plots.

To reproduce the functional test through 1,000 cells (three repetitions each):

```sh
BENCHMARK_CELLS='1 10 100 1000' BENCHMARK_ZORN_MODE=installed \
  BENCHMARK_RSS_INTERVAL=0.02 \
  BENCHMARK_OUT=/husky/henriksson/atrandi/v9_performance_shard/benchmarks/test_up_to_1k \
  Rscript run_performance_shard_benchmark.R
```
