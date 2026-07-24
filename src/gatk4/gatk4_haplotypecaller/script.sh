#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# Stage the reference trio and BAM/BAI into a temp dir using matching
# basenames so GATK can find them
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
staged_bam=$(stage_bam_bai "$tmp_dir" "$par_input" "$par_bai")

# Determine available memory for the JVM
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# unset "false" flags
[[ "$par_disable_optimizations" == "false" ]] && unset par_disable_optimizations
[[ "$par_dont_use_soft_clipped_bases" == "false" ]] && unset par_dont_use_soft_clipped_bases
[[ "$par_floor_blocks" == "false" ]] && unset par_floor_blocks
[[ "$par_force_active" == "false" ]] && unset par_force_active
[[ "$par_sites_only_vcf_output" == "false" ]] && unset par_sites_only_vcf_output
[[ "$par_disable_tool_default_read_filters" == "false" ]] && unset par_disable_tool_default_read_filters
[[ "$par_disable_sequence_dictionary_validation" == "false" ]] && unset par_disable_sequence_dictionary_validation

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args
split_multiple_to_flags "$par_annotation" "--annotation" annotation_args
split_multiple_to_flags "$par_annotation_group" "--annotation-group" annotation_group_args
split_multiple_to_flags "$par_annotations_to_exclude" "--annotations-to-exclude" annotations_to_exclude_args
split_multiple_to_flags "$par_founder_id" "--founder-id" founder_id_args
split_multiple_to_flags "$par_gvcf_gq_bands" "--gvcf-gq-bands" gvcf_gq_bands_args
split_multiple_to_flags "$par_kmer_size" "--kmer-size" kmer_size_args

# Build command arguments array
cmd_args=(
  --input "$staged_bam"
  --output "$par_output"
  --reference "$staged_reference"
  --native-pair-hmm-threads "${meta_cpus:-1}"
  ${par_active_probability_threshold:+--active-probability-threshold "$par_active_probability_threshold"}
  ${par_alleles:+--alleles "$par_alleles"}
  "${annotation_args[@]}"
  "${annotation_group_args[@]}"
  "${annotations_to_exclude_args[@]}"
  ${par_assembly_region_padding:+--assembly-region-padding "$par_assembly_region_padding"}
  ${par_base_quality_score_threshold:+--base-quality-score-threshold "$par_base_quality_score_threshold"}
  ${par_contamination_fraction_to_filter:+--contamination-fraction-to-filter "$par_contamination_fraction_to_filter"}
  ${par_dbsnp:+--dbsnp "$par_dbsnp"}
  ${par_disable_optimizations:+--disable-optimizations}
  ${par_dont_use_soft_clipped_bases:+--dont-use-soft-clipped-bases}
  ${par_emit_ref_confidence:+--emit-ref-confidence "$par_emit_ref_confidence"}
  ${par_floor_blocks:+--floor-blocks}
  ${par_force_active:+--force-active}
  "${founder_id_args[@]}"
  "${gvcf_gq_bands_args[@]}"
  ${par_heterozygosity:+--heterozygosity "$par_heterozygosity"}
  ${par_heterozygosity_stdev:+--heterozygosity-stdev "$par_heterozygosity_stdev"}
  ${par_indel_heterozygosity:+--indel-heterozygosity "$par_indel_heterozygosity"}
  ${par_intervals:+--intervals "$par_intervals"}
  "${kmer_size_args[@]}"
  ${par_max_alternate_alleles:+--max-alternate-alleles "$par_max_alternate_alleles"}
  ${par_max_assembly_region_size:+--max-assembly-region-size "$par_max_assembly_region_size"}
  ${par_max_reads_per_alignment_start:+--max-reads-per-alignment-start "$par_max_reads_per_alignment_start"}
  ${par_min_assembly_region_size:+--min-assembly-region-size "$par_min_assembly_region_size"}
  ${par_min_base_quality_score:+--min-base-quality-score "$par_min_base_quality_score"}
  ${par_min_pruning:+--min-pruning "$par_min_pruning"}
  ${par_minimum_mapping_quality:+--minimum-mapping-quality "$par_minimum_mapping_quality"}
  ${par_output_mode:+--output-mode "$par_output_mode"}
  ${par_pcr_indel_model:+--pcr-indel-model "$par_pcr_indel_model"}
  ${par_pedigree:+--pedigree "$par_pedigree"}
  ${par_ploidy_regions:+--ploidy-regions "$par_ploidy_regions"}
  ${par_sample_name:+--sample-name "$par_sample_name"}
  ${par_sample_ploidy:+--sample-ploidy "$par_sample_ploidy"}
  ${par_stand_call_conf:+--standard-min-confidence-threshold-for-calling "$par_stand_call_conf"}
  --tmp-dir "$tmp_dir"
  ${par_interval_padding:+--interval-padding "$par_interval_padding"}
  ${par_sites_only_vcf_output:+--sites-only-vcf-output}
  --create-output-variant-index "$par_create_output_variant_index"
  --create-output-bam-index "$par_create_output_bam_index"
  "${read_filter_args[@]}"
  "${disable_read_filter_args[@]}"
  ${par_disable_tool_default_read_filters:+--disable-tool-default-read-filters}
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
)

# Run GATK HaplotypeCaller
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" HaplotypeCaller "${cmd_args[@]}"
