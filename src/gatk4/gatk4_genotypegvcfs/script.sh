#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
unset_if_false=(
  par_lenient
  par_annotate_with_num_discovered_alleles
  par_call_genotypes
  par_disable_sequence_dictionary_validation
  par_disable_tool_default_annotations
  par_enable_all_annotations
  par_genomicsdb_shared_posixfs_optimizations
  par_genomicsdb_use_bcf_codec
  par_include_non_variant_sites
  par_input_is_somatic
  par_keep_combined_raw_annotations
  par_merge_input_intervals
  par_only_output_calls_starting_in_intervals
  par_sites_only_vcf_output
  par_use_posteriors_to_calculate_qual
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

if [[ -n "$par_dbsnp_index" && -z "$par_dbsnp" ]]; then
  echo "Error: --dbsnp_index requires --dbsnp." >&2
  exit 1
fi
if [[ -n "$par_population_callset_index" && -z "$par_population_callset" ]]; then
  echo "Error: --population_callset_index requires --population_callset." >&2
  exit 1
fi
if [[ -d "$par_variant" && -n "$par_variant_index" ]]; then
  echo "Error: --variant_index cannot be used when --variant is a GenomicsDB workspace directory." >&2
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

# Stage the inputs into a temp dir so GATK can find the companion files
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_genotypegvcfs.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict") || exit 1

# GATK requires a literal `gendb://` prefix to recognize a GenomicsDB
# workspace. The `--variant` argument may be either a single GVCF file or a
# GenomicsDB workspace directory, so detect which one was provided.
if [ -d "$par_variant" ]; then
  variant_arg="gendb://$par_variant"
else
  # GATK reads a plain --variant VCF without an index, unless it must traverse
  # by intervals (any interval argument)
  variant_index_optional=$([[ "$(gatk_uses_intervals)" == "true" ]] && echo false || echo true)
  variant_arg=$(stage_vcf_with_index "$tmp_dir" "$par_variant" "$par_variant_index" variant "$variant_index_optional") || exit 1
fi
if [[ -n "$par_dbsnp" ]]; then
  staged_dbsnp=$(stage_vcf_with_index "$tmp_dir" "$par_dbsnp" "$par_dbsnp_index" dbsnp) || exit 1
fi
if [[ -n "$par_population_callset" ]]; then
  # GATK reads a plain --population_callset VCF without an index
  staged_population_callset=$(stage_vcf_with_index "$tmp_dir" "$par_population_callset" "$par_population_callset_index" population_callset true) || exit 1
fi
stage_interval_arguments "$tmp_dir" interval_args

# Convert semicolon-separated multi-value arguments into repeated flags
split_multiple_to_flags "$par_annotation" "--annotation" annotation_args
split_multiple_to_flags "$par_annotation_group" "--annotation-group" annotation_group_args
split_multiple_to_flags "$par_annotations_to_exclude" "--annotations-to-exclude" annotations_to_exclude_args
split_multiple_to_flags "$par_force_output_intervals" "--force-output-intervals" force_output_intervals_args
split_multiple_to_flags "$par_founder_id" "--founder-id" founder_id_args
split_multiple_to_flags "$par_keep_specific_combined_raw_annotation" "--keep-specific-combined-raw-annotation" keep_specific_combined_raw_annotation_args

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Build command arguments array
cmd_args=(
  ${par_allele_fraction_error:+--allele-fraction-error "$par_allele_fraction_error"}
  ${par_annotate_with_num_discovered_alleles:+--annotate-with-num-discovered-alleles}
  "${annotation_args[@]}"
  "${annotation_group_args[@]}"
  "${annotations_to_exclude_args[@]}"
  ${par_call_genotypes:+--call-genotypes}
  ${staged_dbsnp:+--dbsnp "$staged_dbsnp"}
  ${par_disable_tool_default_annotations:+--disable-tool-default-annotations}
  ${par_enable_all_annotations:+--enable-all-annotations}
  "${force_output_intervals_args[@]}"
  "${founder_id_args[@]}"
  ${par_genomicsdb_max_alternate_alleles:+--genomicsdb-max-alternate-alleles "$par_genomicsdb_max_alternate_alleles"}
  ${par_genomicsdb_shared_posixfs_optimizations:+--genomicsdb-shared-posixfs-optimizations}
  ${par_genomicsdb_use_bcf_codec:+--genomicsdb-use-bcf-codec}
  ${par_genotype_assignment_method:+--genotype-assignment-method "$par_genotype_assignment_method"}
  ${par_heterozygosity:+--heterozygosity "$par_heterozygosity"}
  ${par_heterozygosity_stdev:+--heterozygosity-stdev "$par_heterozygosity_stdev"}
  ${par_include_non_variant_sites:+--include-non-variant-sites}
  ${par_indel_heterozygosity:+--indel-heterozygosity "$par_indel_heterozygosity"}
  ${par_input_is_somatic:+--input-is-somatic}
  "${interval_args[@]}"
  ${par_keep_combined_raw_annotations:+--keep-combined-raw-annotations}
  "${keep_specific_combined_raw_annotation_args[@]}"
  ${par_max_alternate_alleles:+--max-alternate-alleles "$par_max_alternate_alleles"}
  ${par_max_genotype_count:+--max-genotype-count "$par_max_genotype_count"}
  ${par_max_variants_per_shard:+--max-variants-per-shard "$par_max_variants_per_shard"}
  ${par_merge_input_intervals:+--merge-input-intervals}
  ${par_num_reference_samples_if_no_call:+--num-reference-samples-if-no-call "$par_num_reference_samples_if_no_call"}
  ${par_only_output_calls_starting_in_intervals:+--only-output-calls-starting-in-intervals}
  ${par_pedigree:+--pedigree "$par_pedigree"}
  ${staged_population_callset:+--population-callset "$staged_population_callset"}
  ${par_sample_ploidy:+--sample-ploidy "$par_sample_ploidy"}
  ${par_somatic_quality_threshold:+--somatic-quality-threshold "$par_somatic_quality_threshold"}
  ${par_stand_call_conf:+--standard-min-confidence-threshold-for-calling "$par_stand_call_conf"}
  ${par_tumor_lod_to_emit:+--tumor-lod-to-emit "$par_tumor_lod_to_emit"}
  ${par_use_posteriors_to_calculate_qual:+--use-posteriors-to-calculate-qual}
  ${par_variant_output_filtering:+--variant-output-filtering "$par_variant_output_filtering"}
)

# Run GATK GenotypeGVCFs
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" GenotypeGVCFs \
  --variant "$variant_arg" \
  --reference "$staged_reference" \
  --output "$par_output" \
  "${cmd_args[@]}" \
  --tmp-dir "$tmp_dir" \
  ${par_sites_only_vcf_output:+--sites-only-vcf-output} \
  ${par_create_output_variant_index:+--create-output-variant-index "$par_create_output_variant_index"} \
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation} \
  ${par_lenient:+--lenient}

# Move the output index to --output_index, if requested
if [[ -n "$move_output_index" && ! "$gatk_output_index" -ef "$par_output_index" ]]; then
  mv "$gatk_output_index" "$par_output_index"
fi
