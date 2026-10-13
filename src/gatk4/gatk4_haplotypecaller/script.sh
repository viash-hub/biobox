#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

unset_if_false \
  par_disable_optimizations \
  par_disable_tool_default_annotations \
  par_do_not_run_physical_phasing \
  par_floor_blocks \
  par_use_posteriors_to_calculate_qual

# Check the arguments
require_output_index "$par_output_index" "$par_create_output_variant_index" \
  --output_index --create_output_variant_index
require_output_if_created "$par_create_bam_output" "$par_bam_output" \
  --create_bam_output --bam_output
if is_true "$par_create_bam_output"; then
  if [[ "$par_bam_output" == *.cram ]]; then
    echo "Error: --bam_output must be a BAM file, GATK cannot write CRAM here." >&2
    exit 1
  fi
  require_output_index "$par_bam_output_index" "$par_create_output_bam_index" \
    --bam_output_index --create_output_bam_index
fi

setup_tmp_dir

# Stage the inputs in a temp dir so that GATK can find the companion files.
# The read index is passed with --read-index, so it does not need to be staged.
staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict") || exit 1
gatk_read_args "$par_input" "$par_read_index" read_args
vcf_args=()
stage_optional_vcf "$tmp_dir" "$par_alleles" "$par_alleles_index" --alleles vcf_args
stage_optional_vcf "$tmp_dir" "$par_dbsnp" "$par_dbsnp_index" --dbsnp vcf_args
stage_optional_vcf "$tmp_dir" "$par_population_callset" "$par_population_callset_index" --population_callset vcf_args
stage_interval_arguments "$tmp_dir" interval_args

# Convert semicolon-separated arguments to repeated flags
gatk_annotation_args annotation_args
split_multiple_to_flags "$par_annotations_to_exclude" "--annotations-to-exclude" annotations_to_exclude_args
split_multiple_to_flags "$par_founder_id" "--founder-id" founder_id_args
split_multiple_to_flags "$par_gvcf_gq_bands" "--gvcf-gq-bands" gvcf_gq_bands_args
gatk_assembly_args assembly_args
gatk_engine_args "$tmp_dir" engine_args
gatk_vcf_output_args vcf_output_args
gatk_bam_output_args bam_output_args
gatk_read_filter_args read_filter_args

# Build command arguments array
cmd_args=(
  "${read_args[@]}"
  --output "$par_output"
  --reference "$staged_reference"
  --native-pair-hmm-threads "${meta_cpus:-1}"
  "${vcf_args[@]}"
  "${assembly_args[@]}"
  ${par_activeregion_alt_multiplier:+--activeregion-alt-multiplier "$par_activeregion_alt_multiplier"}
  "${annotation_args[@]}"
  "${annotations_to_exclude_args[@]}"
  ${par_contamination_fraction_per_sample_file:+--contamination-fraction-per-sample-file "$par_contamination_fraction_per_sample_file"}
  ${par_contamination_fraction_to_filter:+--contamination-fraction-to-filter "$par_contamination_fraction_to_filter"}
  ${par_denovo_depth_threshold:+--denovo-depth-threshold "$par_denovo_depth_threshold"}
  ${par_denovo_parent_gq_threshold:+--denovo-parent-gq-threshold "$par_denovo_parent_gq_threshold"}
  ${par_disable_optimizations:+--disable-optimizations}
  ${par_disable_tool_default_annotations:+--disable-tool-default-annotations}
  ${par_do_not_run_physical_phasing:+--do-not-run-physical-phasing}
  ${par_emit_ref_confidence:+--emit-ref-confidence "$par_emit_ref_confidence"}
  ${par_floor_blocks:+--floor-blocks}
  "${founder_id_args[@]}"
  "${gvcf_gq_bands_args[@]}"
  ${par_heterozygosity:+--heterozygosity "$par_heterozygosity"}
  ${par_heterozygosity_stdev:+--heterozygosity-stdev "$par_heterozygosity_stdev"}
  ${par_indel_heterozygosity:+--indel-heterozygosity "$par_indel_heterozygosity"}
  "${interval_args[@]}"
  ${par_mapping_quality_threshold_for_genotyping:+--mapping-quality-threshold-for-genotyping "$par_mapping_quality_threshold_for_genotyping"}
  ${par_max_alternate_alleles:+--max-alternate-alleles "$par_max_alternate_alleles"}
  ${par_max_genotype_count:+--max-genotype-count "$par_max_genotype_count"}
  ${par_max_mnp_distance:+--max-mnp-distance "$par_max_mnp_distance"}
  ${par_maximum_mapping_quality:+--maximum-mapping-quality "$par_maximum_mapping_quality"}
  ${par_minimum_mapping_quality:+--minimum-mapping-quality "$par_minimum_mapping_quality"}
  ${par_num_reference_samples_if_no_call:+--num-reference-samples-if-no-call "$par_num_reference_samples_if_no_call"}
  ${par_output_mode:+--output-mode "$par_output_mode"}
  ${par_pedigree:+--pedigree "$par_pedigree"}
  ${par_ploidy_regions:+--ploidy-regions "$par_ploidy_regions"}
  ${par_sample_name:+--sample-name "$par_sample_name"}
  ${par_sample_ploidy:+--sample-ploidy "$par_sample_ploidy"}
  ${par_standard_min_confidence_threshold_for_calling:+--standard-min-confidence-threshold-for-calling "$par_standard_min_confidence_threshold_for_calling"}
  ${par_use_posteriors_to_calculate_qual:+--use-posteriors-to-calculate-qual}
  "${engine_args[@]}"
  "${vcf_output_args[@]}"
  "${bam_output_args[@]}"
  "${read_filter_args[@]}"
)

# Run GATK HaplotypeCaller
run_gatk HaplotypeCaller "${cmd_args[@]}"

# Move the output index to --output_index
if [[ "$par_create_output_variant_index" != "false" ]]; then
  move_output_index "$(gatk_output_vcf_index_path "$par_output")" "$par_output_index"
fi

# Move the BAM output index to --bam_output_index
if is_true "$par_create_bam_output" && [[ "$par_create_output_bam_index" != "false" ]]; then
  move_output_index "$(gatk_output_bam_index_path "$par_bam_output")" "$par_bam_output_index"
fi
