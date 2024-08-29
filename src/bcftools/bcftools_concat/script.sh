#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Unset parameters
unset_if_false=(
    par_allow_overlaps
    par_compact_PS
    par_ligate
    par_ligate_force
    par_ligate_warn
    par_no_version
    par_naive
    par_naive_force
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Execute bcftools concat with the provided arguments
bcftools concat \
    ${par_allow_overlaps:+-a} \
    ${par_compact_PS:+-c} \
    ${par_remove_duplicates:+-d "$par_remove_duplicates"} \
    ${par_file_list:+-f "$par_file_list"} \
    ${par_ligate:+-l} \
    ${par_ligate_force:+--ligate-force} \
    ${par_ligate_warn:+--ligate-warn} \
    ${par_no_version:+--no-version} \
    ${par_naive:+-n} \
    ${par_naive_force:+--naive-force} \
    ${par_output_type:+--O "$par_output_type"} \
    ${par_min_PQ:+-q "$par_min_PQ"} \
    ${par_regions:+-r "$par_regions"} \
    ${par_regions_file:+-R "$par_regions_file"} \
    ${par_regions_overlap:+--regions-overlap "$par_regions_overlap"} \
    ${par_threads:+--threads "$par_threads"} \
    ${par_verbose:+-v "$par_verbose"} \
    -o $par_output \
    $par_input