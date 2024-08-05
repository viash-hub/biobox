#!/bin/bash

## VIASH START
## VIASH END

# Unset parameters
[[ "$par_strand" == "false" ]] && unset par_strand
[[ "$par_specific_strand" == "false" ]] && unset par_specific_strand
[[ "$par_bed" == "false" ]] && unset par_bed
[[ "$par_header" == "false" ]] && unset par_header
[[ "$par_no_buffer" == "false" ]] && unset par_no_buffer



# Execute bedtools merge with the provided arguments
bedtools merge \
    ${par_strand:+-s} \
    ${par_specific_strand:+-S} \
    ${par_bed:+-bed} \
    ${par_header:+-header} \
    ${par_no_buffer:+-nobuf} \
    ${par_distance:+-d "$par_distance"} \
    ${par_columns:+-c "$par_columns"} \
    ${par_operation:+-0 "$par_operation"} \
    ${par_delimeter:+-delim "$par_delimeter"} \
    ${par_precision:+-prec "$par_precision"} \
    -i "$par_input" \
    > "$par_output"
