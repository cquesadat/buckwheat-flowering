#!/usr/bin/env bash
set -euo pipefail

# Stage the five experiment summaries and designs from Pipeline A into B/data.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
A_EXPERIMENTS="${ROOT}/../A_Image_Segmentation_Modelling/Experiments"
B_DATA="${ROOT}/data"

[[ -d "${A_EXPERIMENTS}" ]] || { echo "Missing Pipeline A experiments: ${A_EXPERIMENTS}" >&2; exit 1; }
[[ -f "${B_DATA}/sowing_dates.csv" ]] || { echo "Missing ${B_DATA}/sowing_dates.csv" >&2; exit 1; }
[[ -f "${B_DATA}/manual_exclusions.csv" ]] || { echo "Missing ${B_DATA}/manual_exclusions.csv" >&2; exit 1; }

experiments=(ETH23 ETH24 LfL23 LfL24 ETH24L)
# Validate all sources before copying any of them.
for exp in "${experiments[@]}"; do
  for filename in "Results_Flowering_Perc_${exp}.csv" "Design_${exp}.csv"; do
    source_file="${A_EXPERIMENTS}/${exp}/${filename}"
    [[ -s "${source_file}" ]] || { echo "Missing or empty: ${source_file}" >&2; exit 1; }
  done
done

mkdir -p "${B_DATA}"
for exp in "${experiments[@]}"; do
  cp -- "${A_EXPERIMENTS}/${exp}/Results_Flowering_Perc_${exp}.csv" "${B_DATA}/"
  cp -- "${A_EXPERIMENTS}/${exp}/Design_${exp}.csv" "${B_DATA}/"
  Rscript "${ROOT}/0.1_Add_DAS.R" "${B_DATA}/Results_Flowering_Perc_${exp}.csv"
done
printf 'Staged summaries and designs for %s experiments in %s\n' "${#experiments[@]}" "${B_DATA}"
