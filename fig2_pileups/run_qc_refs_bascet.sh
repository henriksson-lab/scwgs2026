#!/usr/bin/env bash
#SBATCH -A hpc2n2026-186
#SBATCH -c 10
#SBATCH --mem=16G
#SBATCH -t 0-72:00:00
#SBATCH -J qc
#SBATCH --array=0-1

module load GCC/14.3.0 OpenMPI/5.0.8 HDF5/1.14.6

glob_csv() {
  local pattern=$1
  local matches=()
  local restore_nullglob

  if shopt -q nullglob; then
    restore_nullglob='shopt -s nullglob'
  else
    restore_nullglob='shopt -u nullglob'
  fi
  shopt -s nullglob
  matches=( $pattern )
  eval "$restore_nullglob"

  if ((${#matches[@]} == 0)); then
    printf 'No files matched pattern: %s\n' "$pattern" >&2
    exit 1
  fi

  local IFS=,
  printf '%s' "${matches[*]}"
}

INPUTS=(
  "/home/m/mahogny/mystore/atrandi/v9_wgs_mock4_mda/aligned_cell.*.bam"
  "/home/m/mahogny/mystore/atrandi/v9_wgs_mock5_pta/aligned_cell.*.bam"
)

OUTPUTS=(
  "./data/mda/qc_composition.tsv"
  "./data/pta/qc_composition.tsv"
)

~/github/bascet/target/release/bascet qc ref-composition \
  --in="$(glob_csv "${INPUTS[$SLURM_ARRAY_TASK_ID]}")" \
  --out=${OUTPUTS[$SLURM_ARRAY_TASK_ID]} \
  -@ 10
