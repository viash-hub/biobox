#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

unset_if_false \
  par_allow_old_rms_mapping_quality_annotation_data \
  par_call_genotypes \
  par_convert_to_base_pair_resolution \
  par_disable_tool_default_annotations \
  par_drop_somatic_filtering_annotations \
  par_enable_all_annotations \
  par_ignore_variants_starting_outside_interval \
  par_input_is_somatic

# Check the arguments
require_output_index "$par_output_index" "$par_create_output_variant_index" \
  --output_index --create_output_variant_index

# Stage the inputs in a temp dir so that GATK can find the companion files
setup_tmp_dir

staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict") || exit 1
# GATK reads a plain --variant VCF without an index, unless it must traverse
# by intervals (any interval argument)
variant_index_optional=$([[ "$(gatk_uses_intervals)" == "true" ]] && echo false || echo true)
stage_vcfs_with_indexes "$tmp_dir" "$par_variant" "$par_variant_index" "--variant" "--variant" variant_args "$variant_index_optional"
# GATK reads a plain --dbsnp VCF without an index
dbsnp_args=()
stage_optional_vcf "$tmp_dir" "$par_dbsnp" "$par_dbsnp_index" --dbsnp dbsnp_args true
stage_interval_arguments "$tmp_dir" interval_args

# Convert semicolon-separated multi-value arguments into repeated flags
gatk_annotation_args annotation_args
split_multiple_to_flags "$par_annotations_to_exclude" "--annotations-to-exclude" annotations_to_exclude_args
split_multiple_to_flags "$par_founder_id" "--founder-id" founder_id_args

gatk_engine_args "$tmp_dir" engine_args
gatk_vcf_output_args vcf_output_args

# Build command arguments array
cmd_args=(
  "${variant_args[@]}"
  --reference "$staged_reference"
  --output "$par_output"
  ${par_allow_old_rms_mapping_quality_annotation_data:+--allow-old-rms-mapping-quality-annotation-data}
  "${annotation_args[@]}"
  "${annotations_to_exclude_args[@]}"
  ${par_break_bands_at_multiples_of:+--break-bands-at-multiples-of "$par_break_bands_at_multiples_of"}
  ${par_call_genotypes:+--call-genotypes}
  ${par_combine_variants_distance:+--combine-variants-distance "$par_combine_variants_distance"}
  ${par_convert_to_base_pair_resolution:+--convert-to-base-pair-resolution}
  "${dbsnp_args[@]}"
  ${par_denovo_depth_threshold:+--denovo-depth-threshold "$par_denovo_depth_threshold"}
  ${par_denovo_parent_gq_threshold:+--denovo-parent-gq-threshold "$par_denovo_parent_gq_threshold"}
  ${par_disable_tool_default_annotations:+--disable-tool-default-annotations}
  ${par_drop_somatic_filtering_annotations:+--drop-somatic-filtering-annotations}
  ${par_enable_all_annotations:+--enable-all-annotations}
  "${founder_id_args[@]}"
  ${par_ignore_variants_starting_outside_interval:+--ignore-variants-starting-outside-interval}
  ${par_input_is_somatic:+--input-is-somatic}
  "${interval_args[@]}"
  ${par_max_distance:+--max-distance "$par_max_distance"}
  ${par_pedigree:+--pedigree "$par_pedigree"}
  ${par_ref_padding:+--ref-padding "$par_ref_padding"}
  ${par_variant_output_filtering:+--variant-output-filtering "$par_variant_output_filtering"}
  "${vcf_output_args[@]}"
  "${engine_args[@]}"
)

# Run GATK CombineGVCFs
run_gatk CombineGVCFs "${cmd_args[@]}"

# Move the output index to --output_index
if [[ "$par_create_output_variant_index" != "false" ]]; then
  move_output_index "$(gatk_output_vcf_index_path "$par_output")" "$par_output_index"
fi
