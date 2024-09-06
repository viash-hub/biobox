#!/bin/bash

## VIASH START
## VIASH END

set -e

unset_if_false=(
   par_bam
   par_cram
   par_fast
   par_uncompressed
   par_with_header
   par_header_only
   par_no_header
   par_count
   par_unmap
   par_use_index
   par_fetch_pairs
   par_customized_index
   par_no_PG
   par_write_index
   par_remove_B
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done


samtools view \
    ${par_bam:+-b} \
    ${par_cram:+-C} \
    ${par_fast:+--fast} \
    ${par_uncompressed:+-u} \
    ${par_with_header:+--with-header} \
    ${par_header_only:+-H} \
    ${par_no_header:+--no-header} \
    ${par_count:+-c} \
    ${par_output:+-o "$par_output"} \
    ${par_output_unselected:+-U "$par_output_unselected"} \
    ${par_unmap:+-p "$par_unmap"} \
    ${par_fetch_pairs:+-P "$par_fetch_pairs"} \
    ${par_fai_reference:+-t "$par_fai_reference"} \
    ${par_use_index:+-M "$par_use_index"} \
    ${par_region_file:+--region-file "$par_region_file"} \
    ${par_customized_index:+-X} \
    ${par_target_file:+-L "$par_target_file"} \
    ${par_qname_file:+-N "$par_qname_file"} \
    ${par_read_group:+-r "$par_read_group"} \
    ${par_read_group_file:+-R "$par_read_group_file"} \
    ${par_tag:+-d "$par_tag"} \
    ${par_tag_file:+-D "$par_tag_file"} \
    ${par_min_MQ:+-q "$par_min_MQ"} \
    ${par_library:+-l "$par_library"} \
    ${par_min_qlen:+-m "$par_min_qlen"} \
    ${par_expr:+-e "$par_expr"} \
    ${par_require_flags:+-f "$par_require_flags"} \
    ${par_excl_flags:+-F "$par_excl_flags"} \
    ${par_incl_flags:+--rf "$par_incl_flags"} \
    ${par_excl_all_flags:+-G "$par_excl_all_flags"} \
    ${par_subsample:+--subsample "$par_subsample"} \
    ${par_subsample_seed:+--subsample-seed "$par_subsample_seed"} \
    ${par_add_flags:+--add-flags "$par_add_flags"} \
    ${par_remove_flags:+--remove-flags "$par_remove_flags"} \
    ${par_remove_tag:+-x "$par_remove_tag"} \
    ${par_keep_tag:+--keep-tag "$par_keep_tag"} \
    ${par_remove_B:+-B} \
    ${par_sanitize:+-z "$par_sanitize"} \
    ${par_input_fmt_option:+--input-fmt-option "$par_input_fmt_option"} \
    ${par_output_fmt:+-O "$par_output_fmt"} \
    ${par_output_fmt_option:+--output-fmt-option "$par_output_fmt_option"} \
    ${par_reference:+-T "$par_reference"} \
    ${par_write_index:+--write-index} \
    ${par_no_PG:+--no-PG} \
    "$par_input"

exit 0
