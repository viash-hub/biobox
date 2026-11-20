#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# Check the arguments
require_reference_for_cram "$par_reference" --reference "$par_input"
require_intervals

setup_tmp_dir

# Stage the inputs in a temp dir so that GATK can find the companion files.
# Stage the interval files first, so an index count mismatch fails before the
# (possibly slow) indexing of --variant.
stage_interval_arguments "$tmp_dir" interval_args
gatk_read_args "$par_input" "$par_read_index" read_args
staged_variant=$(stage_vcf_with_index "$tmp_dir" "$par_variant" "$par_variant_index" variant) || exit 1
stage_optional_reference "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict" reference_args

gatk_read_filter_args read_filter_args
gatk_engine_args "$tmp_dir" engine_args

# Build command arguments array
cmd_args=(
  "${read_args[@]}"
  --variant "$staged_variant"
  "${interval_args[@]}"
  --output "$par_output"
  "${reference_args[@]}"
  ${par_max_depth_per_sample:+--max-depth-per-sample "$par_max_depth_per_sample"}
  ${par_maximum_population_allele_frequency:+--maximum-population-allele-frequency "$par_maximum_population_allele_frequency"}
  ${par_maximum_mapping_quality:+--maximum-mapping-quality "$par_maximum_mapping_quality"}
  ${par_minimum_mapping_quality:+--minimum-mapping-quality "$par_minimum_mapping_quality"}
  ${par_minimum_population_allele_frequency:+--minimum-population-allele-frequency "$par_minimum_population_allele_frequency"}
  "${read_filter_args[@]}"
  "${engine_args[@]}"
)

# Run GATK GetPileupSummaries
run_gatk GetPileupSummaries "${cmd_args[@]}"
