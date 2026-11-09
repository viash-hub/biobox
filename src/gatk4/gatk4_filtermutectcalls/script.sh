#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset "false" flags
unset_if_false=(
  par_microbial_mode
  par_mitochondria_mode
  par_sites_only_vcf_output
  par_disable_tool_default_read_filters
  par_disable_sequence_dictionary_validation
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

move_output_index=""
if [[ "$par_create_output_variant_index" != "false" ]]; then
  if [[ -z "$par_output_index" ]]; then
    echo "Error: --output_index is required unless --create_output_variant_index is false." >&2
    exit 1
  fi
  move_output_index="true"
  gatk_output_index=$(gatk_output_vcf_index_path "$par_output")
fi

# Stage the inputs into a temp dir so GATK can find the companion files
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_filtermutectcalls.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

stage_interval_files "$tmp_dir" "$par_intervals" "$par_intervals_index" \
  --intervals --intervals intervals_args
stage_interval_files "$tmp_dir" "$par_exclude_intervals" "$par_exclude_intervals_index" \
  --exclude_intervals --exclude-intervals exclude_intervals_args

staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
staged_variant=$(stage_vcf_with_index "$tmp_dir" "$par_variant" "$par_variant_index" variant)

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_contamination_table" "--contamination-table" contamination_table_args
split_multiple_to_flags "$par_tumor_segmentation" "--tumor-segmentation" tumor_segmentation_args
split_multiple_to_flags "$par_orientation_bias_artifact_priors" "--orientation-bias-artifact-priors" orientation_bias_artifact_priors_args
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args

# Build command arguments array
cmd_args=(
  --variant "$staged_variant"
  --output "$par_output"
  --reference "$staged_reference"
  --stats "$par_stats"
  --filtering-stats "$par_filtering_stats"
  "${contamination_table_args[@]}"
  "${tumor_segmentation_args[@]}"
  "${orientation_bias_artifact_priors_args[@]}"
  "${intervals_args[@]}"
  "${exclude_intervals_args[@]}"
  ${par_contamination_estimate:+--contamination-estimate "$par_contamination_estimate"}
  ${par_distance_on_haplotype:+--distance-on-haplotype "$par_distance_on_haplotype"}
  ${par_f_score_beta:+--f-score-beta "$par_f_score_beta"}
  ${par_false_discovery_rate:+--false-discovery-rate "$par_false_discovery_rate"}
  ${par_initial_threshold:+--initial-threshold "$par_initial_threshold"}
  ${par_interval_set_rule:+--interval-set-rule "$par_interval_set_rule"}
  ${par_log_artifact_prior:+--log-artifact-prior "$par_log_artifact_prior"}
  ${par_log_indel_prior:+--log-indel-prior "$par_log_indel_prior"}
  ${par_log_snv_prior:+--log-snv-prior "$par_log_snv_prior"}
  ${par_long_indel_length:+--long-indel-length "$par_long_indel_length"}
  ${par_max_alt_allele_count:+--max-alt-allele-count "$par_max_alt_allele_count"}
  ${par_max_events_in_haplotype:+--max-events-in-haplotype "$par_max_events_in_haplotype"}
  ${par_max_events_in_region:+--max-events-in-region "$par_max_events_in_region"}
  ${par_max_median_fragment_length_difference:+--max-median-fragment-length-difference "$par_max_median_fragment_length_difference"}
  ${par_max_n_ratio:+--max-n-ratio "$par_max_n_ratio"}
  ${par_microbial_mode:+--microbial-mode}
  ${par_min_allele_fraction:+--min-allele-fraction "$par_min_allele_fraction"}
  ${par_min_median_base_quality:+--min-median-base-quality "$par_min_median_base_quality"}
  ${par_min_median_mapping_quality:+--min-median-mapping-quality "$par_min_median_mapping_quality"}
  ${par_min_median_read_position:+--min-median-read-position "$par_min_median_read_position"}
  ${par_min_reads_per_strand:+--min-reads-per-strand "$par_min_reads_per_strand"}
  ${par_min_slippage_length:+--min-slippage-length "$par_min_slippage_length"}
  ${par_mitochondria_mode:+--mitochondria-mode}
  ${par_normal_p_value_threshold:+--normal-p-value-threshold "$par_normal_p_value_threshold"}
  ${par_pcr_slippage_rate:+--pcr-slippage-rate "$par_pcr_slippage_rate"}
  ${par_threshold_strategy:+--threshold-strategy "$par_threshold_strategy"}
  ${par_unique_alt_read_count:+--unique-alt-read-count "$par_unique_alt_read_count"}
  ${par_variant_output_filtering:+--variant-output-filtering "$par_variant_output_filtering"}
  --tmp-dir "$tmp_dir"
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

# Run GATK FilterMutectCalls
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" FilterMutectCalls "${cmd_args[@]}"

# Move the output index to --output_index, if requested
if [[ -n "$move_output_index" && ! "$gatk_output_index" -ef "$par_output_index" ]]; then
  mv "$gatk_output_index" "$par_output_index"
fi
