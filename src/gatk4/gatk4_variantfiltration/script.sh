#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
unset_if_false=(
  par_lenient
  par_apply_allele_specific_filters
  par_filter_not_in_mask
  par_invalidate_previous_filters
  par_invert_filter_expression
  par_invert_genotype_filter_expression
  par_missing_values_evaluate_as_failing
  par_set_filtered_genotype_to_no_call
  par_sites_only_vcf_output
  par_disable_sequence_dictionary_validation
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Check that the mask index is given together with the mask
if [[ -n "$par_mask_index" && -z "$par_mask" ]]; then
  echo "Error: --mask_index requires --mask." >&2
  exit 1
fi

move_output_index=""
if [[ "$par_create_output_variant_index" != "false" ]]; then
  if [[ -z "$par_output_index" ]]; then
    echo "Error: --output_index is required unless --create_output_variant_index is false." >&2
    exit 1
  fi
  move_output_index="true"
  gatk_output_index=$(gatk_output_vcf_index_path "$par_output")
fi

# Stage the --variant VCF, the --mask file and the interval files with their
# indexes and (if provided) the reference trio into a temp dir so GATK can find
# the companion files
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_variantfiltration.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

# GATK reads a plain --variant VCF without an index, but only when it does not
# traverse by intervals (it then needs random access to the file)
variant_index_optional=$([[ "$(gatk_uses_intervals)" == "true" ]] && echo false || echo true)
staged_variant=$(stage_vcf_with_index "$tmp_dir" "$par_variant" "$par_variant_index" variant "$variant_index_optional") || exit 1

# GATK needs an index for every --mask file, also a plain BED or VCF
staged_mask=""
if [[ -n "$par_mask" ]]; then
  staged_mask=$(stage_vcf_with_index "$tmp_dir" "$par_mask" "$par_mask_index" mask) || exit 1
fi

stage_interval_arguments "$tmp_dir" interval_args

reference_args=()
if [[ -n "$par_reference" ]]; then
  staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict") || exit 1
  reference_args=(--reference "$staged_reference")
fi

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Convert semicolon-separated arguments into repeated flags
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
  ${staged_mask:+--mask "$staged_mask"}
  ${par_mask_description:+--mask-description "$par_mask_description"}
  ${par_mask_extension:+--mask-extension "$par_mask_extension"}
  ${par_mask_name:+--mask-name "$par_mask_name"}
  ${par_missing_values_evaluate_as_failing:+--missing-values-evaluate-as-failing}
  ${par_set_filtered_genotype_to_no_call:+--set-filtered-genotype-to-no-call}
)

# Run gatk VariantFiltration
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" VariantFiltration \
  --variant "$staged_variant" \
  --output "$par_output" \
  "${reference_args[@]}" \
  "${filter_expression_args[@]}" \
  "${filter_name_args[@]}" \
  "${genotype_filter_expression_args[@]}" \
  "${genotype_filter_name_args[@]}" \
  "${extra_args[@]}" \
  "${interval_args[@]}" \
  --tmp-dir "$tmp_dir" \
  ${par_sites_only_vcf_output:+--sites-only-vcf-output} \
  ${par_create_output_variant_index:+--create-output-variant-index "$par_create_output_variant_index"} \
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation} \
  ${par_lenient:+--lenient}

# Move the output index to --output_index, if requested
if [[ -n "$move_output_index" && ! "$gatk_output_index" -ef "$par_output_index" ]]; then
  mv "$gatk_output_index" "$par_output_index"
fi
