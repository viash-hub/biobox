#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset "false" flags
unset_if_false=(
  par_create_bam_output
  par_create_f1r2_tar_gz
  par_dont_use_soft_clipped_bases
  par_force_active
  par_force_call_filtered_alleles
  par_genotype_germline_sites
  par_genotype_pon_sites
  par_ignore_itr_artifacts
  par_mitochondria_mode
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

if [[ -n "$par_create_f1r2_tar_gz" && -z "$par_f1r2_tar_gz" ]]; then
  echo "Error: --create_f1r2_tar_gz requires --f1r2_tar_gz." >&2
  exit 1
fi

move_bam_output_index=""
if [[ -n "$par_create_bam_output" ]]; then
  if [[ -z "$par_bam_output" ]]; then
    echo "Error: --create_bam_output requires --bam_output." >&2
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

# Stage the inputs into a temp dir so GATK can find the companion files
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_mutect2.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

stage_bams_bais "$tmp_dir" "$par_input" "$par_bai" staged_bams
stage_interval_files "$tmp_dir" "$par_intervals" "$par_intervals_index" \
  --intervals --intervals intervals_args
stage_interval_files "$tmp_dir" "$par_exclude_intervals" "$par_exclude_intervals_index" \
  --exclude_intervals --exclude-intervals exclude_intervals_args

staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")

input_args=()
for staged_bam in "${staged_bams[@]}"; do
  input_args+=(--input "$staged_bam")
done

resource_args=()
if [[ -n "$par_germline_resource" ]]; then
  staged_germline_resource=$(stage_vcf_with_index "$tmp_dir" "$par_germline_resource" "$par_germline_resource_index" germline_resource)
  resource_args+=(--germline-resource "$staged_germline_resource")
fi
if [[ -n "$par_panel_of_normals" ]]; then
  staged_panel_of_normals=$(stage_vcf_with_index "$tmp_dir" "$par_panel_of_normals" "$par_panel_of_normals_index" panel_of_normals)
  resource_args+=(--panel-of-normals "$staged_panel_of_normals")
fi
if [[ -n "$par_alleles" ]]; then
  staged_alleles=$(stage_vcf_with_index "$tmp_dir" "$par_alleles" "$par_alleles_index" alleles)
  resource_args+=(--alleles "$staged_alleles")
fi

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_normal_sample" "--normal-sample" normal_sample_args
split_multiple_to_flags "$par_annotation" "--annotation" annotation_args
split_multiple_to_flags "$par_annotation_group" "--annotation-group" annotation_group_args
split_multiple_to_flags "$par_annotations_to_exclude" "--annotations-to-exclude" annotations_to_exclude_args
split_multiple_to_flags "$par_kmer_size" "--kmer-size" kmer_size_args
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args

# Build command arguments array
cmd_args=(
  "${input_args[@]}"
  --output "$par_output"
  --reference "$staged_reference"
  --native-pair-hmm-threads "${meta_cpus:-1}"
  "${normal_sample_args[@]}"
  "${resource_args[@]}"
  "${intervals_args[@]}"
  "${exclude_intervals_args[@]}"
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
  ${par_dont_use_soft_clipped_bases:+--dont-use-soft-clipped-bases}
  ${par_downsampling_stride:+--downsampling-stride "$par_downsampling_stride"}
  ${par_f1r2_max_depth:+--f1r2-max-depth "$par_f1r2_max_depth"}
  ${par_f1r2_median_mq:+--f1r2-median-mq "$par_f1r2_median_mq"}
  ${par_f1r2_min_bq:+--f1r2-min-bq "$par_f1r2_min_bq"}
  ${par_force_active:+--force-active}
  ${par_force_call_filtered_alleles:+--force-call-filtered-alleles}
  ${par_genotype_germline_sites:+--genotype-germline-sites}
  ${par_genotype_pon_sites:+--genotype-pon-sites}
  ${par_ignore_itr_artifacts:+--ignore-itr-artifacts}
  ${par_initial_tumor_lod:+--initial-tumor-lod "$par_initial_tumor_lod"}
  ${par_interval_set_rule:+--interval-set-rule "$par_interval_set_rule"}
  "${kmer_size_args[@]}"
  ${par_max_assembly_region_size:+--max-assembly-region-size "$par_max_assembly_region_size"}
  ${par_max_mnp_distance:+--max-mnp-distance "$par_max_mnp_distance"}
  ${par_max_population_af:+--max-population-af "$par_max_population_af"}
  ${par_max_reads_per_alignment_start:+--max-reads-per-alignment-start "$par_max_reads_per_alignment_start"}
  ${par_max_suspicious_reads_per_alignment_start:+--max-suspicious-reads-per-alignment-start "$par_max_suspicious_reads_per_alignment_start"}
  ${par_min_assembly_region_size:+--min-assembly-region-size "$par_min_assembly_region_size"}
  ${par_min_base_quality_score:+--min-base-quality-score "$par_min_base_quality_score"}
  ${par_min_pruning:+--min-pruning "$par_min_pruning"}
  ${par_minimum_allele_fraction:+--minimum-allele-fraction "$par_minimum_allele_fraction"}
  ${par_minimum_mapping_quality:+--minimum-mapping-quality "$par_minimum_mapping_quality"}
  ${par_mitochondria_mode:+--mitochondria-mode}
  ${par_normal_lod:+--normal-lod "$par_normal_lod"}
  ${par_pcr_indel_model:+--pcr-indel-model "$par_pcr_indel_model"}
  ${par_pcr_indel_qual:+--pcr-indel-qual "$par_pcr_indel_qual"}
  ${par_pcr_snv_qual:+--pcr-snv-qual "$par_pcr_snv_qual"}
  ${par_tumor_lod_to_emit:+--tumor-lod-to-emit "$par_tumor_lod_to_emit"}
  --tmp-dir "$tmp_dir"
  ${par_interval_padding:+--interval-padding "$par_interval_padding"}
  ${par_sites_only_vcf_output:+--sites-only-vcf-output}
  ${par_create_output_variant_index:+--create-output-variant-index "$par_create_output_variant_index"}
  ${par_create_output_bam_index:+--create-output-bam-index "$par_create_output_bam_index"}
  ${par_output_cram_version:+--output-cram-version "$par_output_cram_version"}
  "${read_filter_args[@]}"
  "${disable_read_filter_args[@]}"
  ${par_disable_tool_default_read_filters:+--disable-tool-default-read-filters}
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
)

# Run GATK Mutect2
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" Mutect2 "${cmd_args[@]}"

# Mutect2 has no argument to set the path of its stats file. It always writes
# it next to the output VCF, so move it to --output_stats.
if [[ ! "${par_output}.stats" -ef "$par_output_stats" ]]; then
  mv "${par_output}.stats" "$par_output_stats"
fi

# Move indexes if requested
if [[ -n "$move_output_index" && ! "$gatk_output_index" -ef "$par_output_index" ]]; then
  mv "$gatk_output_index" "$par_output_index"
fi
if [[ -n "$move_bam_output_index" && ! "$gatk_bam_output_index" -ef "$par_bam_output_index" ]]; then
  mv "$gatk_bam_output_index" "$par_bam_output_index"
fi
