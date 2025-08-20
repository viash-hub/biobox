#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false flags to prevent them from being passed to bedtools
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

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Convert semicolon-separated input_b files to array
IFS=';' read -ra input_b_array <<< "$par_input_b"

# Build input B arguments array
input_b_args=()
for file in "${input_b_array[@]}"; do
    input_b_args+=(-b "$file")
done

# Execute bedtools closest
bedtools closest \
    -a "$par_input_a" \
    "${input_b_args[@]}" \
    ${par_distance:+-d} \
    ${par_distance_mode:+-D "$par_distance_mode"} \
    ${par_ignore_overlaps:+-io} \
    ${par_ignore_upstream:+-iu} \
    ${par_ignore_downstream:+-id} \
    ${par_force_upstream:+-fu} \
    ${par_force_downstream:+-fd} \
    ${par_strand:+-s} \
    ${par_different_strand:+-S} \
    ${par_k_closest:+-k "$par_k_closest"} \
    ${par_tie_mode:+-t "$par_tie_mode"} \
    ${par_different_names:+-N} \
    > "$par_output"
