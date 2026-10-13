#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

unset_if_false \
  par_apply_allele_specific_filters \
  par_filter_not_in_mask \
  par_invalidate_previous_filters \
  par_invert_filter_expression \
  par_invert_genotype_filter_expression \
  par_missing_values_evaluate_as_failing \
  par_set_filtered_genotype_to_no_call

# Check the arguments
require_output_index "$par_output_index" "$par_create_output_variant_index" \
  --output_index --create_output_variant_index

# Stage the inputs in a temp dir so that GATK can find the companion files
setup_tmp_dir

# GATK reads a plain --variant VCF without an index, but only when it does not
# traverse by intervals (it then needs random access to the file)
variant_index_optional=$([[ "$(gatk_uses_intervals)" == "true" ]] && echo false || echo true)
staged_variant=$(stage_vcf_with_index "$tmp_dir" "$par_variant" "$par_variant_index" variant "$variant_index_optional") || exit 1
# GATK needs an index for every --mask file, also a plain BED or VCF
mask_args=()
stage_optional_vcf "$tmp_dir" "$par_mask" "$par_mask_index" --mask mask_args
stage_interval_arguments "$tmp_dir" interval_args
stage_optional_reference "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict" reference_args

# Convert semicolon-separated multi-value arguments into repeated flags
split_multiple_to_flags "$par_filter_expression" "--filter-expression" filter_expression_args
split_multiple_to_flags "$par_filter_name" "--filter-name" filter_name_args
split_multiple_to_flags "$par_genotype_filter_expression" "--genotype-filter-expression" genotype_filter_expression_args
split_multiple_to_flags "$par_genotype_filter_name" "--genotype-filter-name" genotype_filter_name_args

gatk_engine_args "$tmp_dir" engine_args
gatk_vcf_output_args vcf_output_args

# Build command arguments array
cmd_args=(
  --variant "$staged_variant"
  --output "$par_output"
  "${reference_args[@]}"
  "${filter_expression_args[@]}"
  "${filter_name_args[@]}"
  "${genotype_filter_expression_args[@]}"
  "${genotype_filter_name_args[@]}"
  ${par_apply_allele_specific_filters:+--apply-allele-specific-filters}
  ${par_cluster_size:+--cluster-size "$par_cluster_size"}
  ${par_cluster_window_size:+--cluster-window-size "$par_cluster_window_size"}
  ${par_filter_not_in_mask:+--filter-not-in-mask}
  ${par_invalidate_previous_filters:+--invalidate-previous-filters}
  ${par_invert_filter_expression:+--invert-filter-expression}
  ${par_invert_genotype_filter_expression:+--invert-genotype-filter-expression}
  "${mask_args[@]}"
  ${par_mask_description:+--mask-description "$par_mask_description"}
  ${par_mask_extension:+--mask-extension "$par_mask_extension"}
  ${par_mask_name:+--mask-name "$par_mask_name"}
  ${par_missing_values_evaluate_as_failing:+--missing-values-evaluate-as-failing}
  ${par_set_filtered_genotype_to_no_call:+--set-filtered-genotype-to-no-call}
  ${par_variant_output_filtering:+--variant-output-filtering "$par_variant_output_filtering"}
  "${interval_args[@]}"
  "${vcf_output_args[@]}"
  "${engine_args[@]}"
)

# Run GATK VariantFiltration
run_gatk VariantFiltration "${cmd_args[@]}"

# Move the output index to --output_index
if [[ "$par_create_output_variant_index" != "false" ]]; then
  move_output_index "$(gatk_output_vcf_index_path "$par_output")" "$par_output_index"
fi
