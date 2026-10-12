#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
[[ "$par_use_original_qualities" == "false" ]] && unset par_use_original_qualities
[[ "$par_disable_tool_default_read_filters" == "false" ]] && unset par_disable_tool_default_read_filters
[[ "$par_disable_sequence_dictionary_validation" == "false" ]] && unset par_disable_sequence_dictionary_validation
[[ "$par_lenient" == "false" ]] && unset par_lenient

# Stage the reference trio, known-sites VCFs and interval files into a temp
# dir using matching basenames so GATK can find them. The read index is passed
# with --read-index, so it does not need to be staged.
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_baserecalibrator.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
gatk_read_args "$par_input" "$par_read_index" read_args
stage_vcfs_with_indexes "$tmp_dir" "$par_known_sites" "$par_known_sites_index" \
  --known_sites --known-sites known_sites_args
stage_interval_arguments "$tmp_dir" interval_args

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args
split_multiple_to_flags "$par_inverted_read_filter" "--inverted-read-filter" inverted_read_filter_args

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Build command arguments array
cmd_args=(
  ${par_binary_tag_name:+--binary-tag-name "$par_binary_tag_name"}
  ${par_bqsr_baq_gap_open_penalty:+--bqsr-baq-gap-open-penalty "$par_bqsr_baq_gap_open_penalty"}
  ${par_default_base_qualities:+--default-base-qualities "$par_default_base_qualities"}
  ${par_deletions_default_quality:+--deletions-default-quality "$par_deletions_default_quality"}
  ${par_indels_context_size:+--indels-context-size "$par_indels_context_size"}
  ${par_insertions_default_quality:+--insertions-default-quality "$par_insertions_default_quality"}
  "${interval_args[@]}"
  ${par_low_quality_tail:+--low-quality-tail "$par_low_quality_tail"}
  ${par_maximum_cycle_value:+--maximum-cycle-value "$par_maximum_cycle_value"}
  ${par_mismatches_context_size:+--mismatches-context-size "$par_mismatches_context_size"}
  ${par_mismatches_default_quality:+--mismatches-default-quality "$par_mismatches_default_quality"}
  ${par_preserve_qscores_less_than:+--preserve-qscores-less-than "$par_preserve_qscores_less_than"}
  ${par_quantizing_levels:+--quantizing-levels "$par_quantizing_levels"}
  ${par_use_original_qualities:+--use-original-qualities}
  "${read_filter_args[@]}"
  "${disable_read_filter_args[@]}"
  ${par_disable_tool_default_read_filters:+--disable-tool-default-read-filters}
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
  ${par_lenient:+--lenient}
  ${par_read_validation_stringency:+--read-validation-stringency "$par_read_validation_stringency"}
  "${inverted_read_filter_args[@]}"
)

# Run GATK BaseRecalibrator
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" BaseRecalibrator \
  "${read_args[@]}" \
  --reference "$staged_reference" \
  --output "$par_output" \
  "${known_sites_args[@]}" \
  "${cmd_args[@]}" \
  --tmp-dir "$tmp_dir"
