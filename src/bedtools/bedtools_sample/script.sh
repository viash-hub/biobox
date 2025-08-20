#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_output_bed" == "false" ]] && unset par_output_bed
[[ "$par_uncompressed_bam" == "false" ]] && unset par_uncompressed_bam
[[ "$par_include_header" == "false" ]] && unset par_include_header
[[ "$par_no_buffer" == "false" ]] && unset par_no_buffer

# Build command arguments array
cmd_args=(
  -i "$par_input"
  ${par_number:+-n "$par_number"}
  ${par_seed:+-seed "$par_seed"}
  ${par_strand_requirement:+-s "$par_strand_requirement"}
  ${par_output_bed:+-bed}
  ${par_uncompressed_bam:+-ubam}
  ${par_include_header:+-header}
  ${par_no_buffer:+-nobuf}
  ${par_input_buffer:+-iobuf "$par_input_buffer"}
)

# Execute bedtools sample
bedtools sample "${cmd_args[@]}" > "$par_output"
