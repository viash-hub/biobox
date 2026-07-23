#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
unset_if_false=(
  par_apply_jexl_filters_first
  par_call_genotypes
  par_exclude_filtered
  par_exclude_non_variants
  par_ignore_non_ref_in_types
  par_invert_mendelian_violation
  par_invert_select
  par_keep_original_ac
  par_keep_original_dp
  par_mendelian_violation
  par_preserve_alleles
  par_remove_unused_alternates
  par_set_filtered_gt_to_nocall
  par_sites_only_vcf_output
  par_disable_tool_default_read_filters
  par_disable_sequence_dictionary_validation
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Convert semicolon-separated arguments into repeated flags
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args
split_multiple_to_flags "$par_drop_genotype_annotation" "--drop-genotype-annotation" drop_genotype_annotation_args
split_multiple_to_flags "$par_drop_info_annotation" "--drop-info-annotation" drop_info_annotation_args
split_multiple_to_flags "$par_exclude_ids" "--exclude-ids" exclude_ids_args
split_multiple_to_flags "$par_exclude_sample_expressions" "--exclude-sample-expressions" exclude_sample_expressions_args
split_multiple_to_flags "$par_exclude_sample_name" "--exclude-sample-name" exclude_sample_name_args
split_multiple_to_flags "$par_keep_ids" "--keep-ids" keep_ids_args
split_multiple_to_flags "$par_sample_name" "--sample-name" sample_name_args
split_multiple_to_flags "$par_sample_expressions" "--sample-expressions" sample_expressions_args
split_multiple_to_flags "$par_select_expressions" "--select" select_expressions_args
split_multiple_to_flags "$par_select_genotype_expressions" "--select-genotype-expressions" select_genotype_expressions_args
split_multiple_to_flags "$par_select_type_to_include" "--select-type-to-include" select_type_to_include_args
split_multiple_to_flags "$par_select_type_to_exclude" "--select-type-to-exclude" select_type_to_exclude_args

# Build command arguments array
cmd_args=(
  --variant "$par_variant"
  --output "$par_output"
  ${par_apply_jexl_filters_first:+--apply-jexl-filters-first}
  ${par_call_genotypes:+--call-genotypes}
  ${par_concordance:+--concordance "$par_concordance"}
  ${par_discordance:+--discordance "$par_discordance"}
  "${drop_genotype_annotation_args[@]}"
  "${drop_info_annotation_args[@]}"
  ${par_exclude_filtered:+--exclude-filtered}
  "${exclude_ids_args[@]}"
  ${par_exclude_intervals:+--exclude-intervals "$par_exclude_intervals"}
  ${par_exclude_non_variants:+--exclude-non-variants}
  "${exclude_sample_expressions_args[@]}"
  "${exclude_sample_name_args[@]}"
  ${par_ignore_non_ref_in_types:+--ignore-non-ref-in-types}
  ${par_intervals:+--intervals "$par_intervals"}
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
  "${select_expressions_args[@]}"
  "${select_genotype_expressions_args[@]}"
  ${par_select_random_fraction:+--select-random-fraction "$par_select_random_fraction"}
  "${select_type_to_exclude_args[@]}"
  "${select_type_to_include_args[@]}"
  ${par_set_filtered_gt_to_nocall:+--set-filtered-gt-to-nocall}
)

# Shared "GATK Engine Options" arguments (see common_argument_groups.yaml)
cmd_args+=(
  ${par_interval_padding:+--interval-padding "$par_interval_padding"}
  ${par_sites_only_vcf_output:+--sites-only-vcf-output}
  --create-output-variant-index "$par_create_output_variant_index"
  --create-output-bam-index "$par_create_output_bam_index"
  "${read_filter_args[@]}"
  "${disable_read_filter_args[@]}"
  ${par_disable_tool_default_read_filters:+--disable-tool-default-read-filters}
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
)

# Run gatk SelectVariants
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" SelectVariants "${cmd_args[@]}"
