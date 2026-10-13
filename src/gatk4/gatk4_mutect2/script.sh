#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# Unset "false" flags
unset_if_false \
  par_create_bam_output \
  par_create_f1r2_tar_gz \
  par_disable_tool_default_annotations \
  par_dont_use_soft_clipped_bases \
  par_enable_dynamic_read_disqualification_for_genotyping \
  par_force_active \
  par_force_call_filtered_alleles \
  par_genotype_germline_sites \
  par_genotype_pon_sites \
  par_ignore_itr_artifacts \
  par_mitochondria_mode

# Check the arguments
require_output_index "$par_output_index" "$par_create_output_variant_index" \
  --output_index --create_output_variant_index
require_output_if_created "$par_create_f1r2_tar_gz" "$par_f1r2_tar_gz" --create_f1r2_tar_gz --f1r2_tar_gz
require_output_if_created "$par_create_bam_output" "$par_bam_output" --create_bam_output --bam_output
if is_true "$par_create_bam_output"; then
  if [[ "$par_bam_output" == *.cram ]]; then
    echo "Error: --bam_output must be a BAM file, GATK cannot write CRAM here." >&2
    exit 1
  fi
  require_output_index "$par_bam_output_index" "$par_create_output_bam_index" \
    --bam_output_index --create_output_bam_index
fi

setup_tmp_dir

# Stage the inputs in a temp dir so that GATK can find the companion files
gatk_read_args "$par_input" "$par_read_index" read_args
stage_interval_arguments "$tmp_dir" interval_args
staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict") || exit 1
resource_args=()
stage_optional_vcf "$tmp_dir" "$par_germline_resource" "$par_germline_resource_index" --germline_resource resource_args
stage_optional_vcf "$tmp_dir" "$par_panel_of_normals" "$par_panel_of_normals_index" --panel_of_normals resource_args
stage_optional_vcf "$tmp_dir" "$par_alleles" "$par_alleles_index" --alleles resource_args

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_normal_sample" "--normal-sample" normal_sample_args
split_multiple_to_flags "$par_annotation" "--annotation" annotation_args
split_multiple_to_flags "$par_annotation_group" "--annotation-group" annotation_group_args
split_multiple_to_flags "$par_annotations_to_exclude" "--annotations-to-exclude" annotations_to_exclude_args
split_multiple_to_flags "$par_kmer_size" "--kmer-size" kmer_size_args
gatk_read_filter_args read_filter_args
gatk_vcf_output_args vcf_output_args
gatk_bam_output_args bam_output_args
gatk_engine_args "$tmp_dir" engine_args

# Build command arguments array
cmd_args=(
  "${read_args[@]}"
  --output "$par_output"
  --reference "$staged_reference"
  --native-pair-hmm-threads "${meta_cpus:-1}"
  "${normal_sample_args[@]}"
  "${resource_args[@]}"
  "${interval_args[@]}"
  ${par_create_f1r2_tar_gz:+--f1r2-tar-gz "$par_f1r2_tar_gz"}
  ${par_create_bam_output:+--bam-output "$par_bam_output"}
  ${par_active_probability_threshold:+--active-probability-threshold "$par_active_probability_threshold"}
  ${par_af_of_alleles_not_in_resource:+--af-of-alleles-not-in-resource "$par_af_of_alleles_not_in_resource"}
  "${annotation_args[@]}"
  "${annotation_group_args[@]}"
  "${annotations_to_exclude_args[@]}"
  ${par_assembly_region_padding:+--assembly-region-padding "$par_assembly_region_padding"}
  ${par_bam_writer_type:+--bam-writer-type "$par_bam_writer_type"}
  ${par_base_quality_score_threshold:+--base-quality-score-threshold "$par_base_quality_score_threshold"}
  ${par_callable_depth:+--callable-depth "$par_callable_depth"}
  ${par_disable_tool_default_annotations:+--disable-tool-default-annotations}
  ${par_dont_use_soft_clipped_bases:+--dont-use-soft-clipped-bases}
  ${par_downsampling_stride:+--downsampling-stride "$par_downsampling_stride"}
  ${par_enable_dynamic_read_disqualification_for_genotyping:+--enable-dynamic-read-disqualification-for-genotyping}
  ${par_f1r2_max_depth:+--f1r2-max-depth "$par_f1r2_max_depth"}
  ${par_f1r2_median_mq:+--f1r2-median-mq "$par_f1r2_median_mq"}
  ${par_f1r2_min_bq:+--f1r2-min-bq "$par_f1r2_min_bq"}
  ${par_force_active:+--force-active}
  ${par_force_call_filtered_alleles:+--force-call-filtered-alleles}
  ${par_genotype_germline_sites:+--genotype-germline-sites}
  ${par_genotype_germline_sites_fraction:+--genotype-germline-sites-fraction "$par_genotype_germline_sites_fraction"}
  ${par_genotype_pon_sites:+--genotype-pon-sites}
  ${par_ignore_itr_artifacts:+--ignore-itr-artifacts}
  ${par_initial_tumor_lod:+--initial-tumor-lod "$par_initial_tumor_lod"}
  "${kmer_size_args[@]}"
  ${par_max_assembly_region_size:+--max-assembly-region-size "$par_max_assembly_region_size"}
  ${par_max_mnp_distance:+--max-mnp-distance "$par_max_mnp_distance"}
  ${par_max_read_length:+--max-read-length "$par_max_read_length"}
  ${par_max_population_af:+--max-population-af "$par_max_population_af"}
  ${par_max_reads_per_alignment_start:+--max-reads-per-alignment-start "$par_max_reads_per_alignment_start"}
  ${par_max_suspicious_reads_per_alignment_start:+--max-suspicious-reads-per-alignment-start "$par_max_suspicious_reads_per_alignment_start"}
  ${par_min_assembly_region_size:+--min-assembly-region-size "$par_min_assembly_region_size"}
  ${par_min_base_quality_score:+--min-base-quality-score "$par_min_base_quality_score"}
  ${par_min_pruning:+--min-pruning "$par_min_pruning"}
  ${par_min_read_length:+--min-read-length "$par_min_read_length"}
  ${par_minimum_allele_fraction:+--minimum-allele-fraction "$par_minimum_allele_fraction"}
  ${par_maximum_mapping_quality:+--maximum-mapping-quality "$par_maximum_mapping_quality"}
  ${par_minimum_mapping_quality:+--minimum-mapping-quality "$par_minimum_mapping_quality"}
  ${par_mitochondria_mode:+--mitochondria-mode}
  ${par_normal_lod:+--normal-lod "$par_normal_lod"}
  ${par_pcr_indel_model:+--pcr-indel-model "$par_pcr_indel_model"}
  ${par_pcr_indel_qual:+--pcr-indel-qual "$par_pcr_indel_qual"}
  ${par_pcr_snv_qual:+--pcr-snv-qual "$par_pcr_snv_qual"}
  ${par_tumor_lod_to_emit:+--tumor-lod-to-emit "$par_tumor_lod_to_emit"}
  "${vcf_output_args[@]}"
  "${bam_output_args[@]}"
  "${read_filter_args[@]}"
  "${engine_args[@]}"
)

# Run GATK Mutect2
run_gatk Mutect2 "${cmd_args[@]}"

# Mutect2 has no argument to set the path of its stats file. It always writes
# it next to the output VCF, so move it to --output_stats.
if [[ ! "${par_output}.stats" -ef "$par_output_stats" ]]; then
  mv "${par_output}.stats" "$par_output_stats"
fi

# Move the output indexes to --output_index and --bam_output_index
if [[ "$par_create_output_variant_index" != "false" ]]; then
  move_output_index "$(gatk_output_vcf_index_path "$par_output")" "$par_output_index"
fi
if is_true "$par_create_bam_output" && [[ "$par_create_output_bam_index" != "false" ]]; then
  move_output_index "$(gatk_output_bam_index_path "$par_bam_output")" "$par_bam_output_index"
fi
