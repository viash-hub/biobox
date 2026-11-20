#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

unset_if_false \
  par_apply_jexl_filters_first \
  par_call_genotypes \
  par_exclude_filtered \
  par_exclude_non_variants \
  par_ignore_non_ref_in_types \
  par_invert_mendelian_violation \
  par_invert_select \
  par_keep_original_ac \
  par_keep_original_dp \
  par_mendelian_violation \
  par_preserve_alleles \
  par_remove_unused_alternates \
  par_set_filtered_gt_to_nocall

# Check the arguments
require_output_index "$par_output_index" "$par_create_output_variant_index" \
  --output_index --create_output_variant_index

# Stage the inputs in a temp dir so that GATK can find the companion files
setup_tmp_dir

# GATK reads a plain --variant VCF without an index, unless it must traverse
# by intervals (any interval argument)
variant_index_optional=$([[ "$(gatk_uses_intervals)" == "true" ]] && echo false || echo true)
staged_variant=$(stage_vcf_with_index "$tmp_dir" "$par_variant" "$par_variant_index" variant "$variant_index_optional") || exit 1
comparison_args=()
stage_optional_vcf "$tmp_dir" "$par_concordance" "$par_concordance_index" --concordance comparison_args
stage_optional_vcf "$tmp_dir" "$par_discordance" "$par_discordance_index" --discordance comparison_args
stage_optional_reference "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict" reference_args
stage_interval_arguments "$tmp_dir" interval_args

# Convert semicolon-separated multi-value arguments into repeated flags
split_multiple_to_flags "$par_drop_genotype_annotation" "--drop-genotype-annotation" drop_genotype_annotation_args
split_multiple_to_flags "$par_drop_info_annotation" "--drop-info-annotation" drop_info_annotation_args
split_multiple_to_flags "$par_exclude_ids" "--exclude-ids" exclude_ids_args
split_multiple_to_flags "$par_exclude_sample_expressions" "--exclude-sample-expressions" exclude_sample_expressions_args
split_multiple_to_flags "$par_exclude_sample_name" "--exclude-sample-name" exclude_sample_name_args
split_multiple_to_flags "$par_keep_ids" "--keep-ids" keep_ids_args
split_multiple_to_flags "$par_sample_name" "--sample-name" sample_name_args
split_multiple_to_flags "$par_sample_expressions" "--sample-expressions" sample_expressions_args
split_multiple_to_flags "$par_select" "--select" select_args
split_multiple_to_flags "$par_select_genotype_expressions" "--select-genotype-expressions" select_genotype_expressions_args
split_multiple_to_flags "$par_select_type_to_include" "--select-type-to-include" select_type_to_include_args
split_multiple_to_flags "$par_select_type_to_exclude" "--select-type-to-exclude" select_type_to_exclude_args

gatk_engine_args "$tmp_dir" engine_args
gatk_vcf_output_args vcf_output_args

# Build command arguments array
cmd_args=(
  --variant "$staged_variant"
  --output "$par_output"
  "${reference_args[@]}"
  ${par_apply_jexl_filters_first:+--apply-jexl-filters-first}
  ${par_call_genotypes:+--call-genotypes}
  "${comparison_args[@]}"
  "${drop_genotype_annotation_args[@]}"
  "${drop_info_annotation_args[@]}"
  ${par_exclude_filtered:+--exclude-filtered}
  "${exclude_ids_args[@]}"
  ${par_exclude_non_variants:+--exclude-non-variants}
  "${exclude_sample_expressions_args[@]}"
  "${exclude_sample_name_args[@]}"
  ${par_ignore_non_ref_in_types:+--ignore-non-ref-in-types}
  "${interval_args[@]}"
  ${par_invert_mendelian_violation:+--invert-mendelian-violation}
  ${par_invert_select:+--invertSelect}
  "${keep_ids_args[@]}"
  ${par_keep_original_ac:+--keep-original-ac}
  ${par_keep_original_dp:+--keep-original-dp}
  ${par_max_filtered_genotypes:+--max-filtered-genotypes "$par_max_filtered_genotypes"}
  ${par_max_fraction_filtered_genotypes:+--max-fraction-filtered-genotypes "$par_max_fraction_filtered_genotypes"}
  ${par_max_indel_size:+--max-indel-size "$par_max_indel_size"}
  ${par_max_nocall_fraction:+--max-nocall-fraction "$par_max_nocall_fraction"}
  ${par_max_nocall_number:+--max-nocall-number "$par_max_nocall_number"}
  ${par_mendelian_violation:+--mendelian-violation}
  ${par_mendelian_violation_qual_threshold:+--mendelian-violation-qual-threshold "$par_mendelian_violation_qual_threshold"}
  ${par_min_filtered_genotypes:+--min-filtered-genotypes "$par_min_filtered_genotypes"}
  ${par_min_fraction_filtered_genotypes:+--min-fraction-filtered-genotypes "$par_min_fraction_filtered_genotypes"}
  ${par_min_indel_size:+--min-indel-size "$par_min_indel_size"}
  ${par_pedigree:+--pedigree "$par_pedigree"}
  ${par_preserve_alleles:+--preserve-alleles}
  ${par_remove_fraction_genotypes:+--remove-fraction-genotypes "$par_remove_fraction_genotypes"}
  ${par_remove_unused_alternates:+--remove-unused-alternates}
  ${par_restrict_alleles_to:+--restrict-alleles-to "$par_restrict_alleles_to"}
  "${sample_expressions_args[@]}"
  "${sample_name_args[@]}"
  "${select_args[@]}"
  "${select_genotype_expressions_args[@]}"
  ${par_select_random_fraction:+--select-random-fraction "$par_select_random_fraction"}
  "${select_type_to_exclude_args[@]}"
  "${select_type_to_include_args[@]}"
  ${par_set_filtered_gt_to_nocall:+--set-filtered-gt-to-nocall}
  ${par_variant_output_filtering:+--variant-output-filtering "$par_variant_output_filtering"}
  "${vcf_output_args[@]}"
  "${engine_args[@]}"
)

# Run GATK SelectVariants
run_gatk SelectVariants "${cmd_args[@]}"

# Move the output index to --output_index
if [[ "$par_create_output_variant_index" != "false" ]]; then
  move_output_index "$(gatk_output_vcf_index_path "$par_output")" "$par_output_index"
fi
