#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
sim_dir="$repo_dir/simulate_mock_reads/atrandi_wgs_sim"
output_root="${PERFORMANCE_SHARD_ROOT:-/husky/henriksson/atrandi/v9_performance_shard}"
seed="${SIM_SEED:-42}"
compression_threads="${SIM_COMPRESSION_THREADS:-8}"
read -r -a cell_counts <<< "${SIM_CELL_COUNTS:-1 10 100 1000 5000 10000 50000 276243}"
# ceil(3,000,000 genomic bases / (150 + 105 genomic bases per pair)).
readonly pairs_per_cell=11765

# The full Illumina 25B lane is an exact total, with approximately 1x per cell.
read_pairs_for_cells() {
    if [[ "$1" == 276243 ]]; then
        echo 3250000000
    else
        echo "$(( $1 * pairs_per_cell ))"
    fi
}

if [[ $# -gt 1 || ( $# -eq 1 && "$1" != --dry-run ) ]]; then
    echo "Usage: $0 [--dry-run]" >&2
    exit 2
fi
[[ "$seed" =~ ^[0-9]+$ ]] || { echo 'SIM_SEED must be a nonnegative integer' >&2; exit 2; }
[[ "$compression_threads" =~ ^[1-9][0-9]{0,2}$ ]] && (( compression_threads <= 128 )) || {
    echo 'SIM_COMPRESSION_THREADS must be between 1 and 128 (workers per mate)' >&2; exit 2;
}
for cells in "${cell_counts[@]}"; do
    [[ "$cells" =~ ^[1-9][0-9]{0,5}$ ]] && (( cells <= 331776 )) || {
        echo "Invalid cell count: $cells (expected 1..331776)" >&2
        exit 2
    }
done

if [[ "${1:-}" == --dry-run ]]; then
    printf 'cells\tread_pairs\tseed\toutput_directory\n'
    for cells in "${cell_counts[@]}"; do
        printf '%s\t%s\t%s\t%s/cells_%s\n' "$cells" "$(read_pairs_for_cells "$cells")" "$seed" "$output_root" "$cells"
    done
    exit 0
fi

cargo build --release --locked --manifest-path "$sim_dir/Cargo.toml" --target-dir "$sim_dir/target"
simulator="$sim_dir/target/release/atrandi-wgs-sim"
binary_sha256="$(sha256sum "$simulator" | cut -d ' ' -f 1)"
source_sha256="$(cd "$sim_dir"; sha256sum Cargo.toml Cargo.lock src/main.rs atrandi_barcodes.tsv | sha256sum | cut -d ' ' -f 1)"
rust_version="$(rustc --version)"

mkdir -p -- "$output_root"
output_root="$(cd -- "$output_root" && pwd)"
exec 9>"$output_root/.simulation.lock"
flock -n 9 || { echo "Another simulation is running in $output_root" >&2; exit 1; }

work_dir=''
trap 'if [[ -n "$work_dir" ]]; then rm -rf -- "$work_dir"; fi' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

for cells in "${cell_counts[@]}"; do
    pairs="$(read_pairs_for_cells "$cells")"
    mean_pairs_per_cell="$pairs_per_cell"
    if [[ "$cells" == 276243 ]]; then
        mean_pairs_per_cell="$(awk 'BEGIN { printf "%.9f", 3250000000 / 276243 }')"
    fi
    destination="$output_root/cells_$cells"
    work_dir="$(mktemp -d "$output_root/.cells_$cells.tmp.XXXXXX")"
    printf 'parameter\tvalue\n' > "$work_dir/parameters.tsv"
    printf '%s\t%s\n' \
        cells "$cells" \
        read_pairs "$pairs" \
        read_pairs_per_cell "$mean_pairs_per_cell" \
        genome_bases 3000000 \
        coverage 1 \
        genomic_bases_per_pair 255 \
        read_length 150 \
        seed "$seed" \
        compression BGZF \
        compression_threads_per_mate "$compression_threads" \
        simulator_source_sha256 "$source_sha256" \
        simulator_binary_sha256 "$binary_sha256" \
        rustc "$rust_version" >> "$work_dir/parameters.tsv"

    if [[ -e "$destination" ]]; then
        if cmp -s "$work_dir/parameters.tsv" "$destination/parameters.tsv" &&
            [[ -f "$destination/checksums.sha256" ]] &&
            (cd "$destination" && sha256sum --check --status checksums.sha256); then
            echo "Verified existing dataset: $destination"
            rm -rf -- "$work_dir"
            work_dir=''
            continue
        fi
        echo "Existing dataset differs or is incomplete: $destination. Use a new output root or move it aside." >&2
        exit 1
    fi

    echo "Generating $cells cells, $pairs read pairs -> $destination"
    "$simulator" --cells "$cells" --reads "$pairs" --seed "$seed" \
        --compression-threads "$compression_threads" \
        --out-prefix "$work_dir/reads" 2> "$work_dir/simulator.log"
    cat "$work_dir/simulator.log"
    (cd "$work_dir" && sha256sum reads_R1.fastq.gz reads_R2.fastq.gz > checksums.sha256)
    # Publish only complete datasets; failures leave no final dataset directory.
    mv -T -- "$work_dir" "$destination"
    work_dir=''
done
