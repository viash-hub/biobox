#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# Check that an index is only given with its VCF input
if [[ -n "$par_alleles_index" && -z "$par_alleles" ]]; then
  echo "Error: --alleles_index requires --alleles." >&2
  exit 1
fi
if [[ -n "$par_dbsnp_index" && -z "$par_dbsnp" ]]; then
  echo "Error: --dbsnp_index requires --dbsnp." >&2
  exit 1
fi
if [[ -n "$par_population_callset_index" && -z "$par_population_callset" ]]; then
  echo "Error: --population_callset_index requires --population_callset." >&2
  exit 1
fi

# unset "false" flags
unset_if_false=(
  par_lenient
  par_create_bam_output
  par_disable_optimizations
  par_disable_tool_default_annotations
  par_do_not_run_physical_phasing
  par_dont_use_soft_clipped_bases
  par_floor_blocks
  par_force_active
  par_force_call_filtered_alleles
  par_use_posteriors_to_calculate_qual
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

move_bam_output_index=""
if [[ -n "$par_create_bam_output" ]]; then
  if [[ -z "$par_bam_output" ]]; then
    echo "Error: --create_bam_output requires --bam_output." >&2
    exit 1
  fi
  if [[ "$par_bam_output" == *.cram ]]; then
    echo "Error: --bam_output must be a BAM file, GATK cannot write CRAM here." >&2
    exit 1
  fi
  if [[ "$par_create_output_bam_index" != "false" ]]; then
    if [[ -z "$par_bam_output_index" ]]; then
      echo "Error: --create_bam_output requires --bam_output_index unless --create_output_bam_index is false." >&2
      exit 1
    fi
    move_bam_output_index="true"
    gatk_bam_output_index=$(gatk_output_bam_index_path "$par_bam_output")
  fi
fi

# Stage the reference trio, VCF inputs and interval files into a temp dir using
# matching basenames so GATK can find them
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_haplotypecaller.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
gatk_read_args "$par_input" "$par_read_index" read_args

if [[ -n "$par_alleles" ]]; then
  staged_alleles=$(stage_vcf_with_index "$tmp_dir" "$par_alleles" "$par_alleles_index" alleles) || exit 1
fi
if [[ -n "$par_dbsnp" ]]; then
  staged_dbsnp=$(stage_vcf_with_index "$tmp_dir" "$par_dbsnp" "$par_dbsnp_index" dbsnp) || exit 1
fi
if [[ -n "$par_population_callset" ]]; then
  staged_population_callset=$(stage_vcf_with_index "$tmp_dir" "$par_population_callset" "$par_population_callset_index" population_callset) || exit 1
fi
stage_interval_arguments "$tmp_dir" interval_args

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_annotation" "--annotation" annotation_args
split_multiple_to_flags "$par_annotation_group" "--annotation-group" annotation_group_args
split_multiple_to_flags "$par_annotations_to_exclude" "--annotations-to-exclude" annotations_to_exclude_args
split_multiple_to_flags "$par_founder_id" "--founder-id" founder_id_args
split_multiple_to_flags "$par_gvcf_gq_bands" "--gvcf-gq-bands" gvcf_gq_bands_args
split_multiple_to_flags "$par_kmer_size" "--kmer-size" kmer_size_args
gatk_read_filter_args read_filter_args

# Build command arguments array
cmd_args=(
  "${read_args[@]}"
  --output "$par_output"
  --reference "$staged_reference"
  --native-pair-hmm-threads "${meta_cpus:-1}"
  ${par_active_probability_threshold:+--active-probability-threshold "$par_active_probability_threshold"}
  ${par_activeregion_alt_multiplier:+--activeregion-alt-multiplier "$par_activeregion_alt_multiplier"}
  ${staged_alleles:+--alleles "$staged_alleles"}
  "${annotation_args[@]}"
  "${annotation_group_args[@]}"
  "${annotations_to_exclude_args[@]}"
  ${par_assembly_region_padding:+--assembly-region-padding "$par_assembly_region_padding"}
  ${par_create_bam_output:+--bam-output "$par_bam_output"}
  ${par_bam_writer_type:+--bam-writer-type "$par_bam_writer_type"}
  ${par_base_quality_score_threshold:+--base-quality-score-threshold "$par_base_quality_score_threshold"}
  ${par_contamination_fraction_per_sample_file:+--contamination-fraction-per-sample-file "$par_contamination_fraction_per_sample_file"}
  ${par_contamination_fraction_to_filter:+--contamination-fraction-to-filter "$par_contamination_fraction_to_filter"}
  ${staged_dbsnp:+--dbsnp "$staged_dbsnp"}
  ${par_denovo_depth_threshold:+--denovo-depth-threshold "$par_denovo_depth_threshold"}
  ${par_denovo_parent_gq_threshold:+--denovo-parent-gq-threshold "$par_denovo_parent_gq_threshold"}
  ${par_disable_optimizations:+--disable-optimizations}
  ${par_disable_tool_default_annotations:+--disable-tool-default-annotations}
  ${par_do_not_run_physical_phasing:+--do-not-run-physical-phasing}
  ${par_dont_use_soft_clipped_bases:+--dont-use-soft-clipped-bases}
  ${par_emit_ref_confidence:+--emit-ref-confidence "$par_emit_ref_confidence"}
  ${par_floor_blocks:+--floor-blocks}
  ${par_force_active:+--force-active}
  ${par_force_call_filtered_alleles:+--force-call-filtered-alleles}
  "${founder_id_args[@]}"
  "${gvcf_gq_bands_args[@]}"
  ${par_heterozygosity:+--heterozygosity "$par_heterozygosity"}
  ${par_heterozygosity_stdev:+--heterozygosity-stdev "$par_heterozygosity_stdev"}
  ${par_indel_heterozygosity:+--indel-heterozygosity "$par_indel_heterozygosity"}
  "${interval_args[@]}"
  "${kmer_size_args[@]}"
  ${par_mapping_quality_threshold_for_genotyping:+--mapping-quality-threshold-for-genotyping "$par_mapping_quality_threshold_for_genotyping"}
  ${par_max_alternate_alleles:+--max-alternate-alleles "$par_max_alternate_alleles"}
  ${par_max_assembly_region_size:+--max-assembly-region-size "$par_max_assembly_region_size"}
  ${par_max_genotype_count:+--max-genotype-count "$par_max_genotype_count"}
  ${par_max_mnp_distance:+--max-mnp-distance "$par_max_mnp_distance"}
  ${par_max_reads_per_alignment_start:+--max-reads-per-alignment-start "$par_max_reads_per_alignment_start"}
  ${par_maximum_mapping_quality:+--maximum-mapping-quality "$par_maximum_mapping_quality"}
  ${par_min_assembly_region_size:+--min-assembly-region-size "$par_min_assembly_region_size"}
  ${par_min_base_quality_score:+--min-base-quality-score "$par_min_base_quality_score"}
  ${par_min_pruning:+--min-pruning "$par_min_pruning"}
  ${par_minimum_mapping_quality:+--minimum-mapping-quality "$par_minimum_mapping_quality"}
  ${par_num_reference_samples_if_no_call:+--num-reference-samples-if-no-call "$par_num_reference_samples_if_no_call"}
  ${par_output_mode:+--output-mode "$par_output_mode"}
  ${par_pcr_indel_model:+--pcr-indel-model "$par_pcr_indel_model"}
  ${par_pedigree:+--pedigree "$par_pedigree"}
  ${par_ploidy_regions:+--ploidy-regions "$par_ploidy_regions"}
  ${staged_population_callset:+--population-callset "$staged_population_callset"}
  ${par_sample_name:+--sample-name "$par_sample_name"}
  ${par_sample_ploidy:+--sample-ploidy "$par_sample_ploidy"}
  ${par_standard_min_confidence_threshold_for_calling:+--standard-min-confidence-threshold-for-calling "$par_standard_min_confidence_threshold_for_calling"}
  ${par_use_posteriors_to_calculate_qual:+--use-posteriors-to-calculate-qual}
  --tmp-dir "$tmp_dir"
  ${par_sites_only_vcf_output:+--sites-only-vcf-output}
  ${par_create_output_variant_index:+--create-output-variant-index "$par_create_output_variant_index"}
  ${par_create_output_bam_index:+--create-output-bam-index "$par_create_output_bam_index"}
  ${par_output_cram_version:+--output-cram-version "$par_output_cram_version"}
  "${read_filter_args[@]}"
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
  ${par_lenient:+--lenient}
)

# Run GATK HaplotypeCaller
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" HaplotypeCaller "${cmd_args[@]}"

# Move indexes if requested
if [[ -n "$move_output_index" && ! "$gatk_output_index" -ef "$par_output_index" ]]; then
  mv "$gatk_output_index" "$par_output_index"
fi
if [[ -n "$move_bam_output_index" && ! "$gatk_bam_output_index" -ef "$par_bam_output_index" ]]; then
  mv "$gatk_bam_output_index" "$par_bam_output_index"
fi
