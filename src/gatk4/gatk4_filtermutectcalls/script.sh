#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# Unset "false" flags
unset_if_false \
  par_microbial_mode \
  par_mitochondria_mode

# Check the arguments
require_output_index "$par_output_index" "$par_create_output_variant_index" \
  --output_index --create_output_variant_index

setup_tmp_dir

# Stage the inputs in a temp dir so that GATK can find the companion files
stage_interval_arguments "$tmp_dir" interval_args
staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict") || exit 1
# GATK reads a plain --variant VCF without an index, unless it must traverse
# by intervals (any interval argument)
variant_index_optional=$([[ "$(gatk_uses_intervals)" == "true" ]] && echo false || echo true)
staged_variant=$(stage_vcf_with_index "$tmp_dir" "$par_variant" "$par_variant_index" variant "$variant_index_optional") || exit 1

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_contamination_table" "--contamination-table" contamination_table_args
split_multiple_to_flags "$par_tumor_segmentation" "--tumor-segmentation" tumor_segmentation_args
split_multiple_to_flags "$par_orientation_bias_artifact_priors" "--orientation-bias-artifact-priors" orientation_bias_artifact_priors_args
gatk_vcf_output_args vcf_output_args
gatk_engine_args "$tmp_dir" engine_args

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
  "${interval_args[@]}"
  ${par_contamination_estimate:+--contamination-estimate "$par_contamination_estimate"}
  ${par_distance_on_haplotype:+--distance-on-haplotype "$par_distance_on_haplotype"}
  ${par_f_score_beta:+--f-score-beta "$par_f_score_beta"}
  ${par_false_discovery_rate:+--false-discovery-rate "$par_false_discovery_rate"}
  ${par_initial_threshold:+--initial-threshold "$par_initial_threshold"}
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
  "${vcf_output_args[@]}"
  "${engine_args[@]}"
)

# Run GATK FilterMutectCalls
run_gatk FilterMutectCalls "${cmd_args[@]}"

# Move the output index to --output_index
if [[ "$par_create_output_variant_index" != "false" ]]; then
  move_output_index "$(gatk_output_vcf_index_path "$par_output")" "$par_output_index"
fi
