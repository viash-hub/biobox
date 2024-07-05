#!/bin/bash

## VIASH START
## VIASH END

set -e

[[ "$par_remove_dups" == "false" ]] && unset par_remove_dups
[[ "$par_customized_index_file" == "false" ]] && unset par_customized_index_file
[[ "$par_sparse" == "false" ]] && unset par_sparse
[[ "$par_remove_overlaps" == "false" ]] && unset par_remove_overlaps

# change the coverage input from X;X;X to X,X,X
par_coverage=$(echo "$par_coverage" | tr ';' ',')

samtools stats \
    ${par_coverage:+-c "$par_coverage"} \
    ${par_remove_dups:+-d} \
    ${par_required_flag:+-f "$par_required_flag"} \
    ${par_filtering_flag:+-F "$par_filtering_flag"} \
    ${par_GC_depth:+--GC-depth "$par_GC_depth"} \
    ${par_insert_size:+-i "$par_insert_size"} \
    ${par_id:+-I "$par_id"} \
    ${par_read_length:+-l "$par_read_length"} \
    ${par_most_inserts:+-m "$par_most_inserts"} \
    ${par_split_prefix:+-P "$par_split_prefix"} \
    ${par_trim_quality:+-q "$par_trim_quality"} \
    ${par_ref_seq:+-r "$par_ref_seq"} \
    ${par_split:+-S "$par_split"} \
    ${par_target_regions:+-t "$par_target_regions"} \
    ${par_sparse:+-x} \
    ${par_remove_overlaps:+-p} \
    ${par_cov_threshold:+-g "$par_cov_threshold"} \
    ${par_input_fmt_option:+-O "$par_input_fmt_option"} \
    ${par_reference:+-R "$par_reference"} \
    "$par_input" \
    > "$par_output"

exit 0