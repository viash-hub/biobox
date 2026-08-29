#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
unset_if_false=(
  par_avoid_nio
  par_bypass_feature_reader
  par_consolidate
  par_merge_input_intervals
  par_overwrite_existing_genomicsdb_workspace
  par_genomicsdb_shared_posixfs_optimizations
  par_validate_sample_name_map
  par_sites_only_vcf_output
  par_disable_tool_default_read_filters
  par_disable_sequence_dictionary_validation
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Convert semicolon-separated multi-value arguments into repeated flags
split_multiple_to_flags "$par_variant" "--variant" variant_args
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Build command arguments array
# Note: do NOT create the output directory ourselves. GenomicsDBImport
# requires --genomicsdb-workspace-path to be an empty or non-existent
# directory unless --overwrite-existing-genomicsdb-workspace is set, so
# GATK must be the one to create it.
cmd_args=(
  ${par_sample_name_map:+--sample-name-map "$par_sample_name_map"}
  --genomicsdb-workspace-path "$par_genomicsdb_workspace_path"
  --intervals "$par_intervals"
  ${par_avoid_nio:+--avoid-nio true}
  ${par_batch_size:+--batch-size "$par_batch_size"}
  ${par_bypass_feature_reader:+--bypass-feature-reader true}
  ${par_consolidate:+--consolidate true}
  ${par_genomicsdb_segment_size:+--genomicsdb-segment-size "$par_genomicsdb_segment_size"}
  ${par_genomicsdb_shared_posixfs_optimizations:+--genomicsdb-shared-posixfs-optimizations true}
  ${par_genomicsdb_vcf_buffer_size:+--genomicsdb-vcf-buffer-size "$par_genomicsdb_vcf_buffer_size"}
  ${par_header:+--header "$par_header"}
  ${par_max_num_intervals_to_import_in_parallel:+--max-num-intervals-to-import-in-parallel "$par_max_num_intervals_to_import_in_parallel"}
  ${par_merge_contigs_into_num_partitions:+--merge-contigs-into-num-partitions "$par_merge_contigs_into_num_partitions"}
  ${par_merge_input_intervals:+--merge-input-intervals true}
  ${par_overwrite_existing_genomicsdb_workspace:+--overwrite-existing-genomicsdb-workspace true}
  ${par_reader_threads:+--reader-threads "$par_reader_threads"}
  ${par_validate_sample_name_map:+--validate-sample-name-map true}
  ${par_interval_padding:+--interval-padding "$par_interval_padding"}
  ${par_sites_only_vcf_output:+--sites-only-vcf-output}
  --create-output-variant-index "$par_create_output_variant_index"
  --create-output-bam-index "$par_create_output_bam_index"
  "${read_filter_args[@]}"
  "${disable_read_filter_args[@]}"
  ${par_disable_tool_default_read_filters:+--disable-tool-default-read-filters}
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
)

# Run GATK GenomicsDBImport
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" GenomicsDBImport "${variant_args[@]}" "${cmd_args[@]}"
