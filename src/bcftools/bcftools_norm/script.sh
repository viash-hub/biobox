#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Unset parameters
unset_if_false=(
    par_atomize
    par_remove_duplicates
    par_force
    par_no_version
    par_do_not_normalize
    par_strict_filter
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Execute bcftools norm with the provided arguments
bcftools norm \
    ${par_atomize:+--atomize} \
    ${par_atom_overlaps:+--atom-overlaps "$par_atom_overlaps"} \
    ${par_check_ref:+--check-ref "$par_check_ref"} \
    ${par_remove_duplicates:+--remove-duplicates} \
    ${par_rm_dup:+--rm-dup "$par_rm_dup"} \
    ${par_fasta_ref:+--fasta-ref "$par_fasta_ref"} \
    ${par_force:+--force} \
    ${par_keep_sum:+--keep-sum "$par_keep_sum"} \
    ${par_multiallelics:+--multiallelics "$par_multiallelics"} \
    ${par_no_version:+--no-version} \
    ${par_do_not_normalize:+--do-not-normalize} \
    ${par_old_rec_tag:+--old-rec-tag "$par_old_rec_tag"} \
    ${par_regions:+--regions "$par_regions"} \
    ${par_regions_file:+--regions-file "$par_regions_file"} \
    ${par_regions_overlap:+--regions-overlap "$par_regions_overlap"} \
    ${par_site_win:+--site-win "$par_site_win"} \
    ${par_strict_filter:+--strict-filter} \
    ${par_targets:+--targets "$par_targets"} \
    ${par_targets_file:+--targets-file "$par_targets_file"} \
    ${par_targets_overlap:+--targets-overlap "$par_targets_overlap"} \
    ${par_threads:+--threads "$par_threads"} \
    ${par_output_type:+-O "$par_output_type"} \
    -o $par_output \
    $par_input
    
