#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_learnreadorientationmodel.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

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
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" LearnReadOrientationModel "${cmd_args[@]}"
