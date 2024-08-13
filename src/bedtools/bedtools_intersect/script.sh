#!/bin/bash

## VIASH START
## VIASH END

unset_if_false=(
    par_write_a
    par_write_b
    par_left_outer_join
    par_write_overlap
    par_write_overlap_plus
    par_report_A_if_no_overlap
    par_number_of_overlaps_A
    par_report_no_overlaps_A
    par_uncompressed_bam
    par_same_strand
    par_opposite_strand
    par_reciprocal_overlap
    par_either_overlap
    par_split
    par_nonamecheck
    par_sorted
    par_filenames
    par_sortout
    par_bed
    par_no_buffer_output
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done


# Create input array 
IFS=";" read -ra input <<< $par_input_b

bedtools intersect \
    ${par_write_a:+-wa} \
    ${par_write_b:+-wb} \
    ${par_left_outer_join:+-loj} \
    ${par_write_overlap:+-wo} \
    ${par_write_overlap_plus:+-wao} \
    ${par_report_A_if_no_overlap:+-u} \
    ${par_number_of_overlaps_A:+-c} \
    ${par_report_no_overlaps_A:+-v} \
    ${par_uncompressed_bam:+-ubam} \
    ${par_same_strand:+-s} \
    ${par_opposite_strand:+-S} \
    ${par_min_overlap_A:+-f "$par_min_overlap_A"} \
    ${par_min_overlap_B:+-F "$par_min_overlap_B"} \
    ${par_reciprocal_overlap:+-r} \
    ${par_either_overlap:+-e} \
    ${par_split:+-split} \
    ${par_genome:+-g "$par_genome"} \
    ${par_nonamecheck:+-nonamecheck} \
    ${par_sorted:+-sorted} \
    ${par_names:+-names "$par_names"} \
    ${par_filenames:+-filenames} \
    ${par_sortout:+-sortout} \
    ${par_bed:+-bed} \
    ${par_header:+-header} \
    ${par_no_buffer_output:+-nobuf} \
    ${par_io_buffer_size:+-iobuf "$par_io_buffer_size"} \
    -a "$par_input_a" \
    ${par_input_b:+ -b ${input[*]}} \
    > "$par_output"
    