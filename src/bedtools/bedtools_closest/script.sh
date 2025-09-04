#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
unset_if_false=(
  par_distance
  par_ignore_overlaps
  par_ignore_upstream
  par_ignore_downstream
  par_force_upstream
  par_force_downstream
  par_strand
  par_different_strand
  par_different_names
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset "$par"
done

# Convert semicolon-separated input_b files to array
IFS=';' read -ra input_b_array <<< "$par_input_b"

# Build command arguments array
cmd_args=(
  -a "$par_input_a"
  ${par_distance:+-d}
  ${par_distance_mode:+-D "$par_distance_mode"}
  ${par_ignore_overlaps:+-io}
  ${par_ignore_upstream:+-iu}
  ${par_ignore_downstream:+-id}
  ${par_force_upstream:+-fu}
  ${par_force_downstream:+-fd}
  ${par_strand:+-s}
  ${par_different_strand:+-S}
  ${par_k_closest:+-k "$par_k_closest"}
  ${par_tie_mode:+-t "$par_tie_mode"}
  ${par_different_names:+-N}
)

# Add multiple input_b files
for file in "${input_b_array[@]}"; do
  cmd_args+=(-b "$file")
done

# Execute bedtools closest
bedtools closest "${cmd_args[@]}" > "$par_output"
