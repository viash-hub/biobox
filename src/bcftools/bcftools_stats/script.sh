#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Unset parameters
unset_if_false=(
    par_1st_allele_only
    par_split_by_ID
    par_verbose 
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Execute bedtools bamtofastq with the provided arguments
bcftools stats \
    -o "$par_output" \
    ${par_1st_allele_only:+--1st-allele-only} \
    ${par_split_by_ID:+--split-by-ID} \
    ${par_verbose:+--verbose} \
    ${par_allele_frequency_bins:+--af-bins "${par_allele_frequency_bins}"} \
    ${par_allele_frequency_tag:+--af-tag "${par_allele_frequency_tag}"} \
    ${par_collapse:+-c "${par_collapse}"} \
    ${par_depth:+-d "${par_depth}"} \
    ${par_exclude:+-e "${par_exclude}"} \
    ${par_exons:+-E "${par_exons}"} \
    ${par_apply_filters:+-f "${par_apply_filters}"} \
    ${par_fasta_reference:+-F "${par_fasta_reference}"} \
    ${par_include:+-i "${par_include}"} \
    ${par_regions:+-r "${par_regions}"} \
    ${par_regions_file:+-R "${par_regions_file}"} \
    ${par_regions_overlap:+--regions-overlap "${par_regions_overlap}"} \
    ${par_samples:+-s "${par_samples}"} \
    ${par_samples_file:+-S "${par_samples_file}"} \
    ${par_targets:+-t "${par_targets}"} \
    ${par_targets_file:+-T "${par_targets_file}"} \
    ${par_targets_overlaps:+--targets-overlap "${par_targets_overlaps}"} \
    ${par_user_tstv:+-u "${par_user_tstv}"} \
    $par_input \

    