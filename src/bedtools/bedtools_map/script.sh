#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags (using loop for many parameters)
unset_if_false=(
  par_reciprocal
  par_either
  par_same_strand
  par_opposite_strand
  par_split
  par_bed_output
  par_header
  par_no_name_check
  par_no_buffer
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset "$par"
done

# Build command arguments array
cmd_args=(
  -a "$par_input_a"
  -b "$par_input_b"
  ${par_columns:+-c "$par_columns"}
  ${par_operations:+-o "$par_operations"}
  ${par_delimiter:+-delim "$par_delimiter"}
  ${par_precision:+-prec "$par_precision"}
  ${par_min_overlap_a:+-f "$par_min_overlap_a"}
  ${par_min_overlap_b:+-F "$par_min_overlap_b"}
  ${par_reciprocal:+-r}
  ${par_either:+-e}
  ${par_same_strand:+-s}
  ${par_opposite_strand:+-S}
  ${par_split:+-split}
  ${par_bed_output:+-bed}
  ${par_header:+-header}
  ${par_genome:+-g "$par_genome"}
  ${par_no_name_check:+-nonamecheck}
  ${par_no_buffer:+-nobuf}
  ${par_io_buffer:+-iobuf "$par_io_buffer"}
)

# Execute bedtools map
bedtools map "${cmd_args[@]}" > "$par_output"
