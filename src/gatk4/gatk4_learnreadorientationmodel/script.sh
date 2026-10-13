#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

setup_tmp_dir

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_input" "--input" input_args

# Build command arguments array
cmd_args=(
  "${input_args[@]}"
  --output "$par_output"
  ${par_convergence_threshold:+--convergence-threshold "$par_convergence_threshold"}
  ${par_max_depth:+--max-depth "$par_max_depth"}
  ${par_num_em_iterations:+--num-em-iterations "$par_num_em_iterations"}
  --tmp-dir "$tmp_dir"
)

# Run GATK LearnReadOrientationModel
run_gatk LearnReadOrientationModel "${cmd_args[@]}"
