#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false flags to prevent them from being passed to bedtools
unset_if_false=(
    par_strand
    par_bed
    par_header
    par_no_buffer
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Execute bedtools merge
bedtools merge \
    -i "$par_input" \
    ${par_strand:+-s} \
    ${par_specific_strand:+-S "$par_specific_strand"} \
    ${par_distance:+-d "$par_distance"} \
    ${par_columns:+-c "$par_columns"} \
    ${par_operation:+-o "$par_operation"} \
    ${par_delimiter:+-delim "$par_delimiter"} \
    ${par_precision:+-prec "$par_precision"} \
    ${par_bed:+-bed} \
    ${par_header:+-header} \
    ${par_no_buffer:+-nobuf} \
    > "$par_output"
