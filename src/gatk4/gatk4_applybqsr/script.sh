#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
unset_if_false=(
  par_use_original_qualities
  par_allow_missing_read_group
  par_emit_original_quals
  par_round_down_quantized
  par_sites_only_vcf_output
  par_disable_tool_default_read_filters
  par_disable_sequence_dictionary_validation
)
for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Stage the BAM/BAI, and the reference if provided, into a temp dir using
# matching basenames so GATK can find them.
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

staged_bam=$(stage_bam_bai "$tmp_dir" "$par_input" "$par_bai")

reference_args=()
if [[ -n "$par_reference" ]]; then
  staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
  reference_args=(--reference "$staged_reference")
fi

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Convert semicolon-separated static quantized quality levels to repeated flags
split_multiple_to_flags "$par_static_quantized_quals" "--static-quantized-quals" static_quantized_quals_args

# --quantize-quals and --static-quantized-quals are mutually exclusive.
# Since --quantize-quals defaults to 0 (i.e. disabled), only pass it
# through when it is explicitly non-zero so that it does not conflict with
# --static-quantized-quals when both are set.
if [[ -n "$par_quantize_quals" && "$par_quantize_quals" != "0" && -n "$par_static_quantized_quals" ]]; then
  echo "Error: --quantize_quals and --static_quantized_quals are mutually exclusive, please set only one of them." >&2
  exit 1
fi

quantize_quals_args=()
if [[ -n "$par_quantize_quals" && "$par_quantize_quals" != "0" ]]; then
  quantize_quals_args=(--quantize-quals "$par_quantize_quals")
fi

# Convert semicolon-separated read-filter arguments to repeated flags
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args

# Build command arguments array
cmd_args=(
  --input "$staged_bam"
  "${reference_args[@]}"
  --bqsr-recal-file "$par_bqsr_recal_file"
  --output "$par_output"
  ${par_allow_missing_read_group:+--allow-missing-read-group}
  ${par_emit_original_quals:+--emit-original-quals}
  ${par_exclude_intervals:+--exclude-intervals "$par_exclude_intervals"}
  ${par_global_qscore_prior:+--global-qscore-prior "$par_global_qscore_prior"}
  ${par_intervals:+--intervals "$par_intervals"}
  ${par_use_original_qualities:+--use-original-qualities}
  ${par_preserve_qscores_less_than:+--preserve-qscores-less-than "$par_preserve_qscores_less_than"}
  ${par_round_down_quantized:+--round-down-quantized}
  "${static_quantized_quals_args[@]}"
  "${quantize_quals_args[@]}"
  --tmp-dir "$tmp_dir"
  ${par_interval_padding:+--interval-padding "$par_interval_padding"}
  ${par_sites_only_vcf_output:+--sites-only-vcf-output}
  --create-output-variant-index "$par_create_output_variant_index"
  --create-output-bam-index "$par_create_output_bam_index"
  "${read_filter_args[@]}"
  "${disable_read_filter_args[@]}"
  ${par_disable_tool_default_read_filters:+--disable-tool-default-read-filters}
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
)

# Run GATK ApplyBQSR
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" ApplyBQSR "${cmd_args[@]}"
