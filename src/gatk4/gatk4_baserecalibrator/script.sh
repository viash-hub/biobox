#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

unset_if_false par_use_original_qualities

setup_tmp_dir

# Stage the inputs in a temp dir so that GATK can find the companion files.
# The read index is passed with --read-index, so it does not need to be staged.
staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict") || exit 1
gatk_read_args "$par_input" "$par_read_index" read_args
stage_vcfs_with_indexes "$tmp_dir" "$par_known_sites" "$par_known_sites_index" \
  --known_sites --known-sites known_sites_args
stage_interval_arguments "$tmp_dir" interval_args

gatk_engine_args "$tmp_dir" engine_args
gatk_read_filter_args read_filter_args

# Build command arguments array
cmd_args=(
  "${read_args[@]}"
  --reference "$staged_reference"
  --output "$par_output"
  "${known_sites_args[@]}"
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
  "${engine_args[@]}"
  "${read_filter_args[@]}"
)

# Run GATK BaseRecalibrator
run_gatk BaseRecalibrator "${cmd_args[@]}"
