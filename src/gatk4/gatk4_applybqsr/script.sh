#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

unset_if_false \
  par_use_original_qualities \
  par_allow_missing_read_group \
  par_emit_original_quals \
  par_round_down_quantized

# Check the arguments
require_reference_for_cram "$par_reference" --reference "$par_input" "$par_output"
require_output_index "$par_output_index" "$par_create_output_bam_index" \
  --output_index --create_output_bam_index

# --quantize-quals and --static-quantized-quals are mutually exclusive.
# Since --quantize-quals defaults to 0 (i.e. disabled), only pass it
# through when it is explicitly non-zero so that it does not conflict with
# --static-quantized-quals when both are set.
if [[ -n "$par_quantize_quals" && "$par_quantize_quals" != "0" && -n "$par_static_quantized_quals" ]]; then
  echo "Error: --quantize_quals and --static_quantized_quals are mutually exclusive, please set only one of them." >&2
  exit 1
fi

setup_tmp_dir

# Stage the inputs in a temp dir so that GATK can find the companion files.
# The read index is passed with --read-index, so it does not need to be staged.
gatk_read_args "$par_input" "$par_read_index" read_args
stage_interval_arguments "$tmp_dir" interval_args
stage_optional_reference "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict" reference_args

# Convert semicolon-separated static quantized quality levels to repeated flags
split_multiple_to_flags "$par_static_quantized_quals" "--static-quantized-quals" static_quantized_quals_args

quantize_quals_args=()
if [[ -n "$par_quantize_quals" && "$par_quantize_quals" != "0" ]]; then
  quantize_quals_args=(--quantize-quals "$par_quantize_quals")
fi

gatk_engine_args "$tmp_dir" engine_args
gatk_bam_output_args bam_output_args
gatk_read_filter_args read_filter_args

# Build command arguments array
cmd_args=(
  "${read_args[@]}"
  "${reference_args[@]}"
  --bqsr-recal-file "$par_bqsr_recal_file"
  --output "$par_output"
  ${par_allow_missing_read_group:+--allow-missing-read-group}
  ${par_emit_original_quals:+--emit-original-quals}
  ${par_global_qscore_prior:+--global-qscore-prior "$par_global_qscore_prior"}
  "${interval_args[@]}"
  ${par_use_original_qualities:+--use-original-qualities}
  ${par_preserve_qscores_less_than:+--preserve-qscores-less-than "$par_preserve_qscores_less_than"}
  ${par_round_down_quantized:+--round-down-quantized}
  "${static_quantized_quals_args[@]}"
  "${quantize_quals_args[@]}"
  "${engine_args[@]}"
  "${bam_output_args[@]}"
  "${read_filter_args[@]}"
)

# Run GATK ApplyBQSR
run_gatk ApplyBQSR "${cmd_args[@]}"

# Move the output index to --output_index
if [[ "$par_create_output_bam_index" != "false" ]]; then
  move_output_index "$(gatk_output_bam_index_path "$par_output")" "$par_output_index"
fi
