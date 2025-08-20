#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_strand" == "false" ]] && unset par_strand
[[ "$par_bed" == "false" ]] && unset par_bed
[[ "$par_header" == "false" ]] && unset par_header
[[ "$par_no_buffer" == "false" ]] && unset par_no_buffer

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
