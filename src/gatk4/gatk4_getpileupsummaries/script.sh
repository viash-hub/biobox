#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset "false" flags
unset_if_false=(
  par_sites_only_vcf_output
  par_disable_tool_default_read_filters
  par_disable_sequence_dictionary_validation
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Stage the interval files, the BAM/BAI, the --variant VCF with its index and
# (if provided) the reference trio into a temp dir so GATK can find the
# companion files
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_getpileupsummaries.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

# Stage the interval files first, so an index count mismatch fails before the
# (possibly slow) indexing of --variant
stage_interval_files "$tmp_dir" "$par_intervals" "$par_intervals_index" \
  --intervals --intervals intervals_args
stage_interval_files "$tmp_dir" "$par_exclude_intervals" "$par_exclude_intervals_index" \
  --exclude_intervals --exclude-intervals exclude_intervals_args

staged_bam=$(stage_bam_bai "$tmp_dir" "$par_input" "$par_bai")
staged_variant=$(stage_vcf_with_index "$tmp_dir" "$par_variant" "$par_variant_index" variant)

reference_args=()
if [[ -n "$par_reference" ]]; then
  staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
  reference_args=(--reference "$staged_reference")
fi

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Convert semicolon-separated arguments to repeated flags
split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_args
split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_args

# Build command arguments array
cmd_args=(
  --input "$staged_bam"
  --variant "$staged_variant"
  "${intervals_args[@]}"
  --output "$par_output"
  "${reference_args[@]}"
  "${exclude_intervals_args[@]}"
  ${par_interval_set_rule:+--interval-set-rule "$par_interval_set_rule"}
  ${par_max_depth_per_sample:+--max-depth-per-sample "$par_max_depth_per_sample"}
  ${par_maximum_population_allele_frequency:+--maximum-population-allele-frequency "$par_maximum_population_allele_frequency"}
  ${par_minimum_mapping_quality:+--minimum-mapping-quality "$par_minimum_mapping_quality"}
  ${par_minimum_population_allele_frequency:+--minimum-population-allele-frequency "$par_minimum_population_allele_frequency"}
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

# Run GATK GetPileupSummaries
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" GetPileupSummaries "${cmd_args[@]}"
