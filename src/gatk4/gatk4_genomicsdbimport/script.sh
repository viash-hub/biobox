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

if [[ -n "$par_variant_index" && -z "$par_variant" ]]; then
  echo "Error: --variant_index requires --variant." >&2
  exit 1
fi
if [[ -n "$par_header_index" && -z "$par_header" ]]; then
  echo "Error: --header_index requires --header." >&2
  exit 1
fi
if [[ -n "$par_sample_name" && -n "$par_sample_name_map" ]]; then
  echo "Error: --sample_name cannot be used together with --sample_name_map." >&2
  exit 1
fi
if [[ -n "$par_sample_name" && -z "$par_variant" ]]; then
  echo "Error: --sample_name requires --variant." >&2
  exit 1
fi
# GATK fails with "/ by zero" when --header is used with --variant files
# without a sample name map
if [[ -n "$par_header" && -n "$par_variant" && -z "$par_sample_name" ]]; then
  echo "Error: --header with --variant requires --sample_name." >&2
  exit 1
fi

# Stage the inputs into a temp dir so GATK can find the companion files
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_genomicsdbimport.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

stage_vcfs_with_indexes "$tmp_dir" "$par_variant" "$par_variant_index" "--variant" "--variant" variant_args

# Write a sample name map for the staged --variant files and pass it to GATK
# in place of --variant
sample_name_map="$par_sample_name_map"
if [[ -n "$par_sample_name" ]]; then
  IFS=';' read -ra sample_names <<< "$par_sample_name"
  if [[ ${#sample_names[@]} -ne $(( ${#variant_args[@]} / 2 )) ]]; then
    echo "Error: --sample_name must be given once for each file in --variant (got ${#sample_names[@]} sample name(s) and $(( ${#variant_args[@]} / 2 )) file(s))." >&2
    exit 1
  fi
  sample_name_map="$tmp_dir/sample_name_map.tsv"
  for i in "${!sample_names[@]}"; do
    printf '%s\t%s\n' "${sample_names[$i]}" "${variant_args[$(( i * 2 + 1 ))]}"
  done > "$sample_name_map"
  variant_args=()
fi
if [[ -n "$par_header" ]]; then
  staged_header=$(stage_vcf_with_index "$tmp_dir" "$par_header" "$par_header_index" header) || exit 1
fi

# Convert semicolon-separated multi-value arguments into repeated flags
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
  ${sample_name_map:+--sample-name-map "$sample_name_map"}
  --genomicsdb-workspace-path "$par_genomicsdb_workspace_path"
  --intervals "$par_intervals"
  ${par_avoid_nio:+--avoid-nio true}
  ${par_batch_size:+--batch-size "$par_batch_size"}
  ${par_bypass_feature_reader:+--bypass-feature-reader true}
  ${par_consolidate:+--consolidate true}
  ${par_genomicsdb_segment_size:+--genomicsdb-segment-size "$par_genomicsdb_segment_size"}
  ${par_genomicsdb_shared_posixfs_optimizations:+--genomicsdb-shared-posixfs-optimizations true}
  ${par_genomicsdb_vcf_buffer_size:+--genomicsdb-vcf-buffer-size "$par_genomicsdb_vcf_buffer_size"}
  ${staged_header:+--header "$staged_header"}
  ${par_max_num_intervals_to_import_in_parallel:+--max-num-intervals-to-import-in-parallel "$par_max_num_intervals_to_import_in_parallel"}
  ${par_merge_contigs_into_num_partitions:+--merge-contigs-into-num-partitions "$par_merge_contigs_into_num_partitions"}
  ${par_merge_input_intervals:+--merge-input-intervals true}
  ${par_overwrite_existing_genomicsdb_workspace:+--overwrite-existing-genomicsdb-workspace true}
  ${par_reader_threads:+--reader-threads "$par_reader_threads"}
  ${par_validate_sample_name_map:+--validate-sample-name-map true}
  ${par_interval_padding:+--interval-padding "$par_interval_padding"}
  ${par_sites_only_vcf_output:+--sites-only-vcf-output}
  ${par_create_output_variant_index:+--create-output-variant-index "$par_create_output_variant_index"}
  ${par_create_output_bam_index:+--create-output-bam-index "$par_create_output_bam_index"}
  ${par_output_cram_version:+--output-cram-version "$par_output_cram_version"}
  "${read_filter_args[@]}"
  "${disable_read_filter_args[@]}"
  ${par_disable_tool_default_read_filters:+--disable-tool-default-read-filters}
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
)

# Run GATK GenomicsDBImport
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" GenomicsDBImport \
  "${variant_args[@]}" \
  "${cmd_args[@]}" \
  --tmp-dir "$tmp_dir"
