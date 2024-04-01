#!/bin/bash

## VIASH START
## VIASH END

set -e

[[ "$par_uncompressed" == "false" ]] && unset par_uncompressed
[[ "$par_minimiser" == "false" ]] && unset par_minimiser
[[ "$par_not_reverse" == "false" ]] && unset par_not_reverse
[[ "$par_homopolymers" == "false" ]] && unset par_homopolymers
[[ "$par_natural_sort" == "false" ]] && unset par_natural_sort
[[ "$par_ascii_sort" == "false" ]] && unset par_ascii_sort
[[ "$par_write_index" == "false" ]] && unset par_write_index
[[ "$par_no_PG" == "false" ]] && unset par_no_PG
[[ "$par_template_coordinate" == "false" ]] && unset par_template_coordinate

samtools sort \
    ${par_compression:+-l "$par_compression"} \
    ${par_uncompressed:+-u} \
    ${par_minimiser:+-M} \
    ${par_not_reverse:+-R} \
    ${par_kmer_size:+-K "$par_kmer_size"} \
    ${par_order:+-I "$par_order"} \
    ${par_window:+-w "$par_window"} \
    ${par_homopolymers:+-H} \
    ${par_natural_sort:+-n} \
    ${par_ascii_sort:+-N} \
    ${par_tag:+-t "$par_tag"} \
    ${par_input_fmt_option:+--input-fmt-option "$par_input_fmt_option"} \
    ${par_template_coordinate:+--template-coordinate} \
    ${par_write_index:+--write-index} \
    ${par_prefix:+-T "$par_prefix"} \
    ${par_no_PG:+--no-PG} \
    ${par_output_fmt:+-O "$par_output_fmt"} \
    ${par_output_fmt_option:+--output-fmt-option "$par_output_fmt_option"} \
    ${par_reference:+--reference "$par_reference"} \
    -o "$par_output" \
    "$par_input"