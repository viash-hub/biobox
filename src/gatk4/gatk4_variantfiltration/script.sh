#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
unset_if_false=(
  par_apply_allele_specific_filters
  par_filter_not_in_mask
  par_invalidate_previous_filters
  par_invert_filter_expression
  par_invert_genotype_filter_expression
  par_missing_values_evaluate_as_failing
  par_set_filtered_genotype_to_no_call
  par_sites_only_vcf_output
  par_disable_tool_default_read_filters
  par_disable_sequence_dictionary_validation
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Compute available memory for the JVM
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Convert semicolon-separated arguments into repeated flags
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args
split_multiple_to_flags "$par_filter_expression" "--filter-expression" filter_expression_args
split_multiple_to_flags "$par_filter_name" "--filter-name" filter_name_args
split_multiple_to_flags "$par_genotype_filter_expression" "--genotype-filter-expression" genotype_filter_expression_args
split_multiple_to_flags "$par_genotype_filter_name" "--genotype-filter-name" genotype_filter_name_args

# Additional general GATK arguments
extra_args=(
  ${par_apply_allele_specific_filters:+--apply-allele-specific-filters}
  ${par_cluster_size:+--cluster-size "$par_cluster_size"}
  ${par_cluster_window_size:+--cluster-window-size "$par_cluster_window_size"}
  ${par_filter_not_in_mask:+--filter-not-in-mask}
  ${par_invalidate_previous_filters:+--invalidate-previous-filters}
  ${par_invert_filter_expression:+--invert-filter-expression}
  ${par_invert_genotype_filter_expression:+--invert-genotype-filter-expression}
  ${par_mask:+--mask "$par_mask"}
  ${par_mask_description:+--mask-description "$par_mask_description"}
  ${par_mask_extension:+--mask-extension "$par_mask_extension"}
  ${par_mask_name:+--mask-name "$par_mask_name"}
  ${par_missing_values_evaluate_as_failing:+--missing-values-evaluate-as-failing}
  ${par_set_filtered_genotype_to_no_call:+--set-filtered-genotype-to-no-call}
)

# Stage the reference trio into a temp dir with matching basenames so GATK
# can find them
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

reference_args=()
if [[ -n "$par_reference" ]]; then
  staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
  reference_args=(--reference "$staged_reference")
fi

# Run gatk VariantFiltration
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" VariantFiltration \
  --variant "$par_input" \
  --output "$par_output" \
  "${reference_args[@]}" \
  "${filter_expression_args[@]}" \
  "${filter_name_args[@]}" \
  "${genotype_filter_expression_args[@]}" \
  "${genotype_filter_name_args[@]}" \
  "${extra_args[@]}" \
  --tmp-dir "$tmp_dir" \
  ${par_interval_padding:+--interval-padding "$par_interval_padding"} \
  ${par_sites_only_vcf_output:+--sites-only-vcf-output} \
  --create-output-variant-index "$par_create_output_variant_index" \
  --create-output-bam-index "$par_create_output_bam_index" \
  "${read_filter_args[@]}" \
  "${disable_read_filter_args[@]}" \
  ${par_disable_tool_default_read_filters:+--disable-tool-default-read-filters} \
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
