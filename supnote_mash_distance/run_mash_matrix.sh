#!/usr/bin/env bash
set -euo pipefail

input_dir="${1:-sim_mash}"
out_dir="${2:-mash_out}"

threads="${THREADS:-8}"
kmer="${KMER:-21}"
sketch_size="${SKETCH_SIZE:-1000}"
min_copies="${MIN_COPIES:-1}"

mkdir -p "$out_dir"

input_list="$out_dir/input_files.txt"
sample_list="$out_dir/sample_ids.txt"
sketch_dir="$out_dir/sketches"
sketch_list="$out_dir/sketch_files.txt"
sketch_prefix="$out_dir/sim_mash_all"
sketch_file="${sketch_prefix}.msh"
edge_file="$out_dir/mash_distances.tsv"
matrix_file="$out_dir/mash_distance_matrix.tsv"
triangle_file="$out_dir/mash_distance_triangle.phylip"

find "$input_dir" -maxdepth 1 \( -type f -o -type l \) \
  \( -name '*.fastq.gz' -o -name '*.fq.gz' -o -name '*.fastq' -o -name '*.fq' \) \
  | sort > "$input_list"

if [[ ! -s "$input_list" ]]; then
  echo "No FASTQ files found in: $input_dir" >&2
  exit 1
fi

mkdir -p "$sketch_dir"
> "$sample_list"
> "$sketch_list"

echo "Sketching $(wc -l < "$input_list") FASTQ files from $input_dir"
while IFS= read -r fastq; do
  sample="$(basename "$fastq")"
  sample="${sample%.fastq.gz}"
  sample="${sample%.fq.gz}"
  sample="${sample%.fastq}"
  sample="${sample%.fq}"

  sketch_path="$sketch_dir/$sample"
  echo "$sample" >> "$sample_list"
  echo "${sketch_path}.msh" >> "$sketch_list"

  mash sketch \
    -r \
    -I "$sample" \
    -k "$kmer" \
    -s "$sketch_size" \
    -m "$min_copies" \
    -o "$sketch_path" \
    "$fastq"
done < "$input_list"

echo "Combining per-file sketches into $sketch_file"
mash paste -l "$sketch_prefix" "$sketch_list"

echo "Writing Mash edge list to $edge_file"
mash triangle -E -p "$threads" "$sketch_file" > "$edge_file"

echo "Writing Mash relaxed PHYLIP triangle to $triangle_file"
mash triangle -p "$threads" "$sketch_file" > "$triangle_file"

echo "Writing square distance matrix to $matrix_file"
awk '
  BEGIN {
    FS = OFS = "\t"
  }
  FNR == NR {
    n++
    labels[n] = $0
    sub(/^.*\//, "", labels[n])
    sub(/\.fastq\.gz$/, "", labels[n])
    sub(/\.fq\.gz$/, "", labels[n])
    sub(/\.fastq$/, "", labels[n])
    sub(/\.fq$/, "", labels[n])
    idx[$0] = n
    dist[n, n] = 0
    next
  }
  {
    i = idx[$1]
    j = idx[$2]
    if (i && j) {
      dist[i, j] = $3
      dist[j, i] = $3
    }
    next
  }
  END {
    printf "sample"
    for (i = 1; i <= n; i++) {
      printf "%s%s", OFS, labels[i]
    }
    printf ORS

    for (i = 1; i <= n; i++) {
      printf "%s", labels[i]
      for (j = 1; j <= n; j++) {
        value = ((i, j) in dist) ? dist[i, j] : "NA"
        printf "%s%s", OFS, value
      }
      printf ORS
    }
  }
' "$sample_list" "$edge_file" > "$matrix_file"

echo "Done"
echo "Sketch: $sketch_file"
echo "Edge list: $edge_file"
echo "Distance matrix: $matrix_file"
