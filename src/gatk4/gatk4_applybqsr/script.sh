#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
unset_if_false=(
  par_lenient
  par_use_original_qualities
  par_allow_missing_read_group
  par_emit_original_quals
  par_round_down_quantized
  par_disable_tool_default_read_filters
  par_disable_sequence_dictionary_validation
)
for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

IFS=';' read -ra input_files <<< "$par_input"
for read_file in "${input_files[@]}" "$par_output"; do
  if [[ "$read_file" == *.cram && -z "$par_reference" ]]; then
    echo "Error: --reference is required for CRAM input or output." >&2
    exit 1
  fi
done

move_output_index=""
if [[ "$par_create_output_bam_index" != "false" ]]; then
  if [[ -z "$par_output_index" ]]; then
    echo "Error: --output_index is required unless --create_output_bam_index is false." >&2
    exit 1
  fi
  move_output_index="true"
  gatk_output_index=$(gatk_output_bam_index_path "$par_output")
fi

# Stage the interval files and the reference (if provided) into a temp dir
# using matching basenames so GATK can find them. The read index is passed
# with --read-index, so it does not need to be staged.
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_applybqsr.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

gatk_read_args "$par_input" "$par_read_index" read_args
stage_interval_arguments "$tmp_dir" interval_args

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
  --tmp-dir "$tmp_dir"
  ${par_create_output_bam_index:+--create-output-bam-index "$par_create_output_bam_index"}
  ${par_output_cram_version:+--output-cram-version "$par_output_cram_version"}
  "${read_filter_args[@]}"
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
  ${par_lenient:+--lenient}
)

# Run GATK ApplyBQSR
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" ApplyBQSR "${cmd_args[@]}"

# Move the output index to --output_index, if requested
if [[ -n "$move_output_index" && ! "$gatk_output_index" -ef "$par_output_index" ]]; then
  mv "$gatk_output_index" "$par_output_index"
fi
