#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Unset parameters
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

# Execute bedtools merge with the provided arguments
bedtools merge \
    ${par_strand:+-s} \
    ${par_specific_strand:+-S "$par_specific_strand"} \
    ${par_bed:+-bed} \
    ${par_header:+-header} \
    ${par_no_buffer:+-nobuf} \
    ${par_distance:+-d "$par_distance"} \
    ${par_columns:+-c "$par_columns"} \
    ${par_operation:+-o "$par_operation"} \
    ${par_delimiter:+-delim "$par_delimiter"} \
    ${par_precision:+-prec "$par_precision"} \
    -i "$par_input" \
    > "$par_output"
