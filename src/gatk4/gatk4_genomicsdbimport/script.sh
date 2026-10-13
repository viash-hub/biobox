#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

unset_if_false \
  par_avoid_nio \
  par_bypass_feature_reader \
  par_consolidate \
  par_merge_input_intervals \
  par_overwrite_existing_genomicsdb_workspace \
  par_genomicsdb_shared_posixfs_optimizations \
  par_validate_sample_name_map

# Check the arguments
require_intervals
if [[ -n "$par_variant_index" && -z "$par_variant" ]]; then
  echo "Error: --variant_index requires --variant." >&2
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

# Stage the inputs in a temp dir so that GATK can find the companion files
setup_tmp_dir

stage_vcfs_with_indexes "$tmp_dir" "$par_variant" "$par_variant_index" "--variant" "--variant" variant_args
stage_interval_arguments "$tmp_dir" interval_args

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
header_args=()
stage_optional_vcf "$tmp_dir" "$par_header" "$par_header_index" --header header_args

gatk_engine_args "$tmp_dir" engine_args

# Build command arguments array
# Note: do NOT create the output directory ourselves. GenomicsDBImport
# requires --genomicsdb-workspace-path to be an empty or non-existent
# directory unless --overwrite-existing-genomicsdb-workspace is set, so
# GATK must be the one to create it.
cmd_args=(
  "${variant_args[@]}"
  ${sample_name_map:+--sample-name-map "$sample_name_map"}
  --genomicsdb-workspace-path "$par_genomicsdb_workspace_path"
  "${interval_args[@]}"
  ${par_avoid_nio:+--avoid-nio}
  ${par_batch_size:+--batch-size "$par_batch_size"}
  ${par_bypass_feature_reader:+--bypass-feature-reader}
  ${par_consolidate:+--consolidate}
  ${par_genomicsdb_segment_size:+--genomicsdb-segment-size "$par_genomicsdb_segment_size"}
  ${par_genomicsdb_shared_posixfs_optimizations:+--genomicsdb-shared-posixfs-optimizations}
  ${par_genomicsdb_vcf_buffer_size:+--genomicsdb-vcf-buffer-size "$par_genomicsdb_vcf_buffer_size"}
  "${header_args[@]}"
  ${par_max_num_intervals_to_import_in_parallel:+--max-num-intervals-to-import-in-parallel "$par_max_num_intervals_to_import_in_parallel"}
  ${par_merge_contigs_into_num_partitions:+--merge-contigs-into-num-partitions "$par_merge_contigs_into_num_partitions"}
  ${par_merge_input_intervals:+--merge-input-intervals}
  ${par_overwrite_existing_genomicsdb_workspace:+--overwrite-existing-genomicsdb-workspace}
  ${par_reader_threads:+--reader-threads "$par_reader_threads"}
  ${par_validate_sample_name_map:+--validate-sample-name-map}
  "${engine_args[@]}"
)

# Run GATK GenomicsDBImport
run_gatk GenomicsDBImport "${cmd_args[@]}"
