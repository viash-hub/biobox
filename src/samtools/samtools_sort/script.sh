#!/bin/bash

## VIASH START
## VIASH END

set -e

unset_if_false=(
   par_uncompressed
   par_minimiser
   par_not_reverse
   par_homopolymers
   par_natural_sort
   par_ascii_sort
   par_template_coordinate
   par_write_index
   par_no_PG
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done


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

# save text files containing the output of samtools view for later comparison
samtools view "$par_output" -o "$par_output".txt