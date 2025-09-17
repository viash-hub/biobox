#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags (using loop for many parameters)
unset_if_false=(
  par_reciprocal
  par_either_overlap
  par_same_strand
  par_opposite_strand
  par_remove_entire
  par_remove_if_all_overlap
  par_write_original_b
  par_write_overlap_counts
  par_output_bed
  par_include_header
  par_use_split
  par_sorted_input
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
  ${par_min_overlap_a:+-f "$par_min_overlap_a"}
  ${par_min_overlap_b:+-F "$par_min_overlap_b"}
  ${par_reciprocal:+-r}
  ${par_either_overlap:+-e}
  ${par_same_strand:+-s}
  ${par_opposite_strand:+-S}
  ${par_remove_entire:+-A}
  ${par_remove_if_all_overlap:+-N}
  ${par_write_original_b:+-wb}
  ${par_write_overlap_counts:+-wo}
  ${par_output_bed:+-bed}
  ${par_include_header:+-header}
  ${par_use_split:+-split}
  ${par_sorted_input:+-sorted}
  ${par_genome_file:+-g "$par_genome_file"}
  ${par_no_name_check:+-nonamecheck}
  ${par_no_buffer:+-nobuf}
  ${par_input_buffer:+-iobuf "$par_input_buffer"}
)

# Execute bedtools subtract and redirect output to the specified output file
bedtools subtract "${cmd_args[@]}" > "$par_output"
