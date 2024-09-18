#!/bin/bash

## VIASH START
## VIASH END

set -e

test_dir="${metal_executable}/test_data"

unset_if_false=(
   par_paired
   par_in_sam
   par_out_sam
   par_spliced_is_unique
   par_per_gene
   par_per_contig
   par_per_cell
   par_no_sort_output
   par_buffer_whole_contig
   par_ignore_umi
   par_subset
   par_log2stderr
   par_read_length
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

umi_tools dedup \
    --stdin "$par_input" \
    ${par_in_sam:+--in-sam} \
    --stdout "$par_output" \
    ${par_out_sam:+--out-sam} \
    ${par_paired:+--paired} \
    ${par_output_stats:+--output-stats "$par_output_stats"} \
    ${par_extract_umi_method:+--extract-umi-method "$par_extract_umi_method"} \
    ${par_umi_tag:+--umi-tag "$par_umi_tag"} \
    ${par_umi_separator:+--umi-separator "$par_umi_separator"} \
    ${par_umi_tag_split:+--umi-tag-split "$par_umi_tag_split"} \
    ${par_umi_tag_delimiter:+--umi-tag-delimiter "$par_umi_tag_delimiter"} \
    ${par_cell_tag:+--cell-tag "$par_cell_tag"} \
    ${par_cell_tag_split:+--cell-tag-split "$par_cell_tag_split"} \
    ${par_cell_tag_delimiter:+--cell-tag-delimiter "$par_cell_tag_delimiter"} \
    ${par_method:+--method "$par_method"} \
    ${par_edit_distance_threshold:+--edit-distance-threshold "$par_edit_distance_threshold"} \
    ${par_spliced_is_unique:+--spliced-is-unique} \
    ${par_soft_clip_threshold:+--soft-clip-threshold "$par_soft_clip_threshold"} \
    ${par_multimapping_detection_method:+--multimapping-detection-method "$par_multimapping_detection_method"} \
    ${par_read_length:+--read-length} \
    ${par_per_gene:+--per-gene} \
    ${par_gene_tag:+--gene-tag "$par_gene_tag"} \
    ${par_assigned_status_tag:+--assigned-status-tag "$par_assigned_status_tag"} \
    ${par_skip_tags_regex:+--skip-tags-regex "$par_skip_tags_regex"} \
    ${par_per_contig:+--per-contig} \
    ${par_gene_transcript_map:+--gene-transcript-map "$par_gene_transcript_map"} \
    ${par_per_cell:+--per-cell} \
    ${par_mapping_quality:+--mapping-quality "$par_mapping_quality"} \
    ${par_unmapped_reads:+--unmapped-reads "$par_unmapped_reads"} \
    ${par_chimeric_pairs:+--chimeric-pairs "$par_chimeric_pairs"} \
    ${par_unpaired_reads:+--unpaired-reads "$par_unpaired_reads"} \
    ${par_ignore_umi:+--ignore-umi} \
    ${par_subset:+--subset "$par_subset"} \
    ${par_chrom:+--chrom "$par_chrom"} \
    ${par_no_sort_output:+--no-sort-output} \
    ${par_buffer_whole_contig:+--buffer-whole-contig} \
    ${par_log:+-L "$par_log"} \
    ${par_log2stderr:+--log2stderr} \
    ${par_verbose:+-v "$par_verbose"} \
    ${par_error:+-E "$par_error"} \
    ${par_temp_dir:+--temp-dir "$par_temp_dir"} \
    ${par_compresslevel:+--compresslevel "$par_compresslevel"} \
    ${par_timeit:+--timeit "$par_timeit"} \
    ${par_timeit_name:+--timeit-name "$par_timeit_name"} \
    ${par_timeit_header:+--timeit-header "$par_timeit_header"} \
    ${par_random_seed:+--random-seed "$par_random_seed"}

exit 0
