#!/bin/bash

## VIASH START
## VIASH END

set -e

unset_if_false=(
  par_no_suffix
  par_suffix
  par_use_oq
  par_copy_tags
  par_casava
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

if [[ "$meta_name" == "samtools_fasta" ]]; then
  subcommand=fasta
elif [[ "$meta_name" == "samtools_fastq" ]]; then
  subcommand=fastq
else
  echo "Unrecognized component name" && exit 1
fi
samtools "$subcommand" \
    ${par_no_suffix:+-n} \
    ${par_suffix:+-N} \
    ${par_use_oq:+-O} \
    ${par_singleton:+-s "$par_singleton"} \
    ${par_copy_tags:+-t} \
    ${par_copy_tags_list:+-T "$par_copy_tags_list"} \
    ${par_read1:+-1 "$par_read1"} \
    ${par_read2:+-2 "$par_read2"} \
    ${par_output_reads:+-o "$par_output_reads"} \
    ${par_output_reads_both:+-0 "$par_output_reads_both"} \
    ${par_filter_flags:+-f "$par_filter_flags"} \
    ${par_excl_flags:+-F "$par_excl_flags"} \
    ${par_incl_flags:+--rf "$par_incl_flags"} \
    ${par_excl_flags_all:+-G "$par_excl_flags_all"} \
    ${par_aux_tag:+-d "$par_aux_tag"} \
    ${par_aux_tag_file:+-D "$par_aux_tag_file"} \
    ${par_casava:+-i} \
    ${par_compression:+-c "$par_compression"} \
    ${par_index1:+--i1 "$par_index1"} \
    ${par_index2:+--i2 "$par_index2"} \
    ${par_barcode_tag:+--barcode-tag "$par_barcode_tag"} \
    ${par_quality_tag:+--quality-tag "$par_quality_tag"} \
    ${par_index_format:+--index-format "$par_index_format"} \
    "$par_input" \
    > "$par_output"

