#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset boolean flags that are false
[[ "$par_reciprocal" == "false" ]] && unset par_reciprocal
[[ "$par_either_overlap" == "false" ]] && unset par_either_overlap
[[ "$par_same_strand" == "false" ]] && unset par_same_strand
[[ "$par_opposite_strand" == "false" ]] && unset par_opposite_strand
[[ "$par_remove_entire" == "false" ]] && unset par_remove_entire
[[ "$par_remove_if_all_overlap" == "false" ]] && unset par_remove_if_all_overlap
[[ "$par_write_original_b" == "false" ]] && unset par_write_original_b
[[ "$par_write_overlap_counts" == "false" ]] && unset par_write_overlap_counts
[[ "$par_output_bed" == "false" ]] && unset par_output_bed
[[ "$par_include_header" == "false" ]] && unset par_include_header
[[ "$par_use_split" == "false" ]] && unset par_use_split
[[ "$par_sorted_input" == "false" ]] && unset par_sorted_input
[[ "$par_no_name_check" == "false" ]] && unset par_no_name_check
[[ "$par_no_buffer" == "false" ]] && unset par_no_buffer

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
