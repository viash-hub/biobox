#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# Unset "false" flags
unset_if_false par_create_tumor_segmentation

# Check the arguments
require_output_if_created "$par_create_tumor_segmentation" "$par_tumor_segmentation" \
  --create_tumor_segmentation --tumor_segmentation

setup_tmp_dir

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
run_gatk CalculateContamination "${cmd_args[@]}"
