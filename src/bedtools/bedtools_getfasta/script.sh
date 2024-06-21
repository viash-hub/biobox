#!/usr/bin/env bash
set -eo pipefail

unset_if_false=( par_rna par_strandedness par_tab par_bed_out par_name par_name_only par_split par_full_header )

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

bedtools getfasta \
    -fi "$par_input_fasta" \
    -bed "$par_input_bed" \
    ${par_rna:+-rna} \
    ${par_name:+-name} \
    ${par_name_only:+-nameOnly} \
    ${par_tab:+-tab} \
    ${par_bed_out:+-bedOut} \
    ${par_strandedness:+-s} \
    ${par_split:+-split} \
    ${par_full_header:+-fullHeader} > "$par_output"

