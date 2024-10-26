#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Unset parameters
unset_if_false=(
    par_full
    par_inheader
    par_outheader
    par_header
    par_ignorecase
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

bedtools groupby \
    ${par_full:+-full} \
    ${par_inheader:+-inheader} \
    ${par_outheader:+-outheader} \
    ${par_header:+-header} \
    ${par_ignorecase:+-ignorecase} \
    ${par_precision:+-prec "$par_precision"} \
    ${par_delimiter:+-delim "$par_delimiter"} \
    -i "$par_input" \
    -g "$par_groupby" \
    -c "$par_column" \
    ${par_operation:+-o "$par_operation"} \
    > "$par_output"
    