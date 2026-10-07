#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

[[ "$par_create_tumor_segmentation" == "false" ]] && unset par_create_tumor_segmentation

if [[ -n "$par_create_tumor_segmentation" && -z "$par_tumor_segmentation" ]]; then
  echo "Error: --create_tumor_segmentation requires --tumor_segmentation." >&2
  exit 1
fi

tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_calculatecontamination.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Build command arguments array
cmd_args=(
  --input "$par_input"
  --output "$par_output"
  ${par_matched_normal:+--matched-normal "$par_matched_normal"}
  ${par_create_tumor_segmentation:+--tumor-segmentation "$par_tumor_segmentation"}
  ${par_high_coverage_ratio_threshold:+--high-coverage-ratio-threshold "$par_high_coverage_ratio_threshold"}
  ${par_low_coverage_ratio_threshold:+--low-coverage-ratio-threshold "$par_low_coverage_ratio_threshold"}
  --tmp-dir "$tmp_dir"
)

# Run GATK CalculateContamination
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" CalculateContamination "${cmd_args[@]}"
