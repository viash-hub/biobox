#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Unset parameters
unset_if_false=(
    par_force
    par_keep_sites
    par_no_version
    par_single_overlaps
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Execute bcftools annotate with the provided arguments
bcftools annotate \
    ${par_annotations:+-a "$par_annotations"} \
    ${par_columns:+-c "$par_columns"} \
    ${par_columns_file:+-C "$par_columns_file"} \
    ${par_exclude:+-e "$par_exclude"} \
    ${par_force:+--force} \
    ${par_header_line:+-H "$par_header_line"} \
    ${par_header_lines:+-h "$par_header_lines"} \
    ${par_set_id:+-I "$par_set_id"} \
    ${par_include:+-i "$par_include"} \
    ${par_keep_sites:+-k} \
    ${par_merge_logic:+-l "$par_merge_logic"} \
    ${par_mark_sites:+-m "$par_mark_sites"} \
    ${par_min_overlap:+--min-overlap "$par_min_overlap"} \
    ${par_no_version:+--no-version} \
    ${par_samples_file:+-S "$par_samples_file"} \
    ${par_output_type:+-O "$par_output_type"} \
    ${par_pair_logic:+--pair-logic "$par_pair_logic"} \
    ${par_regions:+-r "$par_regions"} \
    ${par_regions_file:+-R "$par_regions_file"} \
    ${par_regions_overlap:+--regions-overlap "$par_regions_overlap"} \
    ${par_rename_annotations:+--rename-annotations "$par_rename_annotations"} \
    ${par_rename_chromosomes:+--rename-chromosomes "$par_rename_chromosomes"} \
    ${par_samples:+-s "$par_samples"} \
    ${par_single_overlaps:+--single-overlaps} \
    ${par_threads:+--threads "$par_threads"} \
    ${par_remove:+-x "$par_remove"} \
    -o $par_output \
    $par_input
    