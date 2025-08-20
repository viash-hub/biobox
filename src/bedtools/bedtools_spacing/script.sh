#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_output_bed" == "false" ]] && unset par_output_bed
[[ "$par_include_header" == "false" ]] && unset par_include_header
[[ "$par_no_buffer" == "false" ]] && unset par_no_buffer

# Build command arguments array
cmd_args=(
  -i "$par_input"
  ${par_output_bed:+-bed}
  ${par_include_header:+-header}
  ${par_no_buffer:+-nobuf}
  ${par_input_buffer:+-iobuf "$par_input_buffer"}
)

# Execute bedtools spacing
bedtools spacing "${cmd_args[@]}" > "$par_output"
