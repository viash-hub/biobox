#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

unset_if_false \
  par_allow_old_rms_mapping_quality_annotation_data \
  par_annotate_with_num_discovered_alleles \
  par_call_genotypes \
  par_disable_tool_default_annotations \
  par_enable_all_annotations \
  par_genomicsdb_shared_posixfs_optimizations \
  par_genomicsdb_use_bcf_codec \
  par_include_non_variant_sites \
  par_input_is_somatic \
  par_keep_combined_raw_annotations \
  par_merge_input_intervals \
  par_use_posteriors_to_calculate_qual

# Check the arguments
if [[ -d "$par_variant" && -n "$par_variant_index" ]]; then
  echo "Error: --variant_index cannot be used when --variant is a GenomicsDB workspace directory." >&2
  exit 1
fi
require_output_index "$par_output_index" "$par_create_output_variant_index" \
  --output_index --create_output_variant_index

# Stage the inputs in a temp dir so that GATK can find the companion files
setup_tmp_dir

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
dbsnp_args=()
stage_optional_vcf "$tmp_dir" "$par_dbsnp" "$par_dbsnp_index" --dbsnp dbsnp_args
# GATK reads a plain --population_callset VCF without an index
population_callset_args=()
stage_optional_vcf "$tmp_dir" "$par_population_callset" "$par_population_callset_index" \
  --population_callset population_callset_args true
stage_interval_arguments "$tmp_dir" interval_args

# Convert semicolon-separated multi-value arguments into repeated flags
gatk_annotation_args annotation_args
split_multiple_to_flags "$par_annotations_to_exclude" "--annotations-to-exclude" annotations_to_exclude_args
# Pass --force_output_intervals files first, then --force_output_interval_strings
split_multiple_to_flags "$par_force_output_intervals" "--force-output-intervals" force_output_intervals_args
split_multiple_to_flags "$par_force_output_interval_strings" "--force-output-intervals" force_output_interval_strings_args
split_multiple_to_flags "$par_founder_id" "--founder-id" founder_id_args
split_multiple_to_flags "$par_keep_specific_combined_raw_annotation" "--keep-specific-combined-raw-annotation" keep_specific_combined_raw_annotation_args

gatk_engine_args "$tmp_dir" engine_args
gatk_vcf_output_args vcf_output_args

# Build command arguments array
cmd_args=(
  --variant "$variant_arg"
  --reference "$staged_reference"
  --output "$par_output"
  ${par_allele_fraction_error:+--allele-fraction-error "$par_allele_fraction_error"}
  ${par_allow_old_rms_mapping_quality_annotation_data:+--allow-old-rms-mapping-quality-annotation-data}
  ${par_annotate_with_num_discovered_alleles:+--annotate-with-num-discovered-alleles}
  "${annotation_args[@]}"
  "${annotations_to_exclude_args[@]}"
  ${par_call_genotypes:+--call-genotypes}
  "${dbsnp_args[@]}"
  ${par_denovo_depth_threshold:+--denovo-depth-threshold "$par_denovo_depth_threshold"}
  ${par_denovo_parent_gq_threshold:+--denovo-parent-gq-threshold "$par_denovo_parent_gq_threshold"}
  ${par_disable_tool_default_annotations:+--disable-tool-default-annotations}
  ${par_enable_all_annotations:+--enable-all-annotations}
  "${force_output_intervals_args[@]}"
  "${force_output_interval_strings_args[@]}"
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
  ${par_merge_input_intervals:+--merge-input-intervals}
  ${par_num_reference_samples_if_no_call:+--num-reference-samples-if-no-call "$par_num_reference_samples_if_no_call"}
  ${par_pedigree:+--pedigree "$par_pedigree"}
  "${population_callset_args[@]}"
  ${par_sample_ploidy:+--sample-ploidy "$par_sample_ploidy"}
  ${par_somatic_quality_threshold:+--somatic-quality-threshold "$par_somatic_quality_threshold"}
  ${par_standard_min_confidence_threshold_for_calling:+--standard-min-confidence-threshold-for-calling "$par_standard_min_confidence_threshold_for_calling"}
  ${par_tumor_lod_to_emit:+--tumor-lod-to-emit "$par_tumor_lod_to_emit"}
  ${par_use_posteriors_to_calculate_qual:+--use-posteriors-to-calculate-qual}
  ${par_variant_output_filtering:+--variant-output-filtering "$par_variant_output_filtering"}
  "${vcf_output_args[@]}"
  "${engine_args[@]}"
)

# Run GATK GenotypeGVCFs
run_gatk GenotypeGVCFs "${cmd_args[@]}"

# Move the output index to --output_index
if [[ "$par_create_output_variant_index" != "false" ]]; then
  move_output_index "$(gatk_output_vcf_index_path "$par_output")" "$par_output_index"
fi
