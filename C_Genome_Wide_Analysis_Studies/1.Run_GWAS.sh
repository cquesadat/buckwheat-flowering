#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
B="${ROOT}/../B_Flowering_dynamics_Modelling/BLUEs_BLUPs"
command -v Rscript >/dev/null || { echo 'Rscript is required' >&2; exit 1; }
[[ -d "$B" ]] || { echo "Missing Pipeline B BLUE/BLUP directory: $B" >&2; exit 1; }
cd "$ROOT"
shopt -s nullglob
for pair in 'ETH23:ETH23' 'ETH24:ETH24' 'ALL:ETH_across_years' 'BBCH_ETH23:BBCH_2023' 'BBCH_ETH24:BBCH_2024' 'BBCH_ALL:BBCH_across_years'; do
  dataset="${pair%%:*}"
  source_group="${pair#*:}"
  source_dir="$B/$source_group"
  [[ -d "$source_dir" ]] || { echo "Missing B phenotype group: $source_dir" >&2; exit 1; }
  files=("$source_dir/"*_BLUEs_BLUPs.csv)
  if [[ "$dataset" == 'BBCH_ALL' ]]; then
    files=("$source_dir/"BLUPs_BLUEs_*.csv)
    if ((${#files[@]} == 0)); then
      files=("$source_dir/"*_BLUEs_BLUPs.csv)
    fi
  fi
  ((${#files[@]} > 0)) || { echo "No phenotype files: $source_dir" >&2; exit 1; }
  for file in "${files[@]}"; do
    name="$(basename "$file")"
    if [[ "$dataset" == 'BBCH_ALL' && "$name" == BLUPs_BLUEs_* ]]; then
      trait="${name#BLUPs_BLUEs_}"
      trait="${trait%.csv}"
    else
      trait="${name%_BLUEs_BLUPs.csv}"
    fi
    echo "GWAS: $dataset / $trait"
    Rscript 1.1.GWAS_batch.R "$dataset" "$trait" 3 "$file"
  done
done
