#!/usr/bin/env bash
set -euo pipefail

# GORG-Dark assembled SAGs -> Bascet TIRP for 31-mer counting.
# 150-base fragments, 30-base overlap. No host filtering is performed here;
# filtered_nohost.1.tirp.gz is the requested downstream input filename.
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
gorg_dir="${GORG_DIR:-/husky/henriksson/atrandi/v9_wgs_gorg}"
gorg_input="${GORG_INPUT:-${gorg_dir}/GORG-Dark_SAGs.zip}"
gorg_output="${gorg_dir}/filtered_nohost.1.tirp.gz"
gorg_url="https://ndownloader.figshare.com/files/54602522"
gorg_md5="780173f01b6cc923c5709978966eb6a5"

for path in "$gorg_output" "${gorg_output}.tbi" "${gorg_output}.records.tsv"; do
    if [[ -e "$path" || -L "$path" ]]; then
        printf 'Refusing to overwrite existing output: %s\n' "$path" >&2
        exit 1
    fi
done

mkdir -p -- "$gorg_dir"
if [[ -z "${GORG_INPUT:-}" ]]; then
    if [[ ! -f "$gorg_input" ]]; then
        printf 'Downloading GORG-Dark SAG assemblies to %s\n' "$gorg_input"
        curl --fail --location --retry 5 --retry-delay 5 \
            --continue-at - --output "${gorg_input}.part" "$gorg_url"
        printf '%s  %s\n' "$gorg_md5" "${gorg_input}.part" | md5sum --check -
        mv -- "${gorg_input}.part" "$gorg_input"
    else
        printf '%s  %s\n' "$gorg_md5" "$gorg_input" | md5sum --check -
    fi
elif [[ ! -e "$gorg_input" ]]; then
    printf 'GORG_INPUT does not exist: %s\n' "$gorg_input" >&2
    exit 1
fi

cargo run --release --locked --manifest-path "${repo_dir}/sag2tirp/Cargo.toml" -- \
    --input "$gorg_input" \
    --output "$gorg_output" \
    --fragment-length 150 \
    --kmer-length 31 \
    --phred 0

bgzip -t "$gorg_output"
printf 'Indexed cells: '
tabix -l "$gorg_output" | wc -l
printf 'Output: %s\n' "$gorg_output"
