#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Unset variables that are false
unset_if_false=( par_histogram par_depth_per_position par_counts_only par_mean_depth 
                 par_same_strand par_different_strand par_reciprocal par_either 
                 par_split par_bed_output par_header par_sorted par_no_name_check )

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Build input B arguments array from semicolon-separated string
input_b_args=()
IFS=';' read -ra input_b_files <<< "$par_input_b"
for file in "${input_b_files[@]}"; do
    input_b_args+=(-b "$file")
done

# Execute bedtools coverage
bedtools coverage \
    -a "$par_input_a" \
    "${input_b_args[@]}" \
    ${par_histogram:+-hist} \
    ${par_depth_per_position:+-d} \
    ${par_counts_only:+-counts} \
    ${par_mean_depth:+-mean} \
    ${par_same_strand:+-s} \
    ${par_different_strand:+-S} \
    ${par_min_overlap_a:+-f "$par_min_overlap_a"} \
    ${par_min_overlap_b:+-F "$par_min_overlap_b"} \
    ${par_reciprocal:+-r} \
    ${par_either:+-e} \
    ${par_split:+-split} \
    ${par_bed_output:+-bed} \
    ${par_header:+-header} \
    ${par_sorted:+-sorted} \
    ${par_genome:+-g "$par_genome"} \
    ${par_no_name_check:+-nonamecheck} \
    > "$par_output"
