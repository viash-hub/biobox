#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
[[ "$par_use_original_qualities" == "false" ]] && unset par_use_original_qualities
[[ "$par_sites_only_vcf_output" == "false" ]] && unset par_sites_only_vcf_output
[[ "$par_disable_tool_default_read_filters" == "false" ]] && unset par_disable_tool_default_read_filters
[[ "$par_disable_sequence_dictionary_validation" == "false" ]] && unset par_disable_sequence_dictionary_validation

# Stage the reference trio and BAM/BAI into a temp dir using matching
# basenames so GATK can find them
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
staged_bam=$(stage_bam_bai "$tmp_dir" "$par_input" "$par_bai")

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_known_sites" "--known-sites" known_sites_args
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args

# Determine available memory for the JVM
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Build command arguments array
cmd_args=(
  ${par_binary_tag_name:+--binary-tag-name "$par_binary_tag_name"}
  ${par_bqsr_baq_gap_open_penalty:+--bqsr-baq-gap-open-penalty "$par_bqsr_baq_gap_open_penalty"}
  ${par_default_base_qualities:+--default-base-qualities "$par_default_base_qualities"}
  ${par_deletions_default_quality:+--deletions-default-quality "$par_deletions_default_quality"}
  ${par_indels_context_size:+--indels-context-size "$par_indels_context_size"}
  ${par_insertions_default_quality:+--insertions-default-quality "$par_insertions_default_quality"}
  ${par_intervals:+--intervals "$par_intervals"}
  ${par_low_quality_tail:+--low-quality-tail "$par_low_quality_tail"}
  ${par_maximum_cycle_value:+--maximum-cycle-value "$par_maximum_cycle_value"}
  ${par_mismatches_context_size:+--mismatches-context-size "$par_mismatches_context_size"}
  ${par_mismatches_default_quality:+--mismatches-default-quality "$par_mismatches_default_quality"}
  ${par_preserve_qscores_less_than:+--preserve-qscores-less-than "$par_preserve_qscores_less_than"}
  ${par_quantizing_levels:+--quantizing-levels "$par_quantizing_levels"}
  ${par_use_original_qualities:+--use-original-qualities}
  ${par_interval_padding:+--interval-padding "$par_interval_padding"}
  ${par_sites_only_vcf_output:+--sites-only-vcf-output}
  --create-output-variant-index "$par_create_output_variant_index"
  --create-output-bam-index "$par_create_output_bam_index"
  "${read_filter_args[@]}"
  "${disable_read_filter_args[@]}"
  ${par_disable_tool_default_read_filters:+--disable-tool-default-read-filters}
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
)

# Run GATK BaseRecalibrator
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" BaseRecalibrator \
  --input "$staged_bam" \
  --reference "$staged_reference" \
  --output "$par_output" \
  "${known_sites_args[@]}" \
  "${cmd_args[@]}" \
  --tmp-dir "$tmp_dir"
