#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
unset_if_false=(
  par_lenient
  par_allow_old_rms_mapping_quality_annotation_data
  par_call_genotypes
  par_convert_to_base_pair_resolution
  par_disable_tool_default_annotations
  par_drop_somatic_filtering_annotations
  par_enable_all_annotations
  par_ignore_variants_starting_outside_interval
  par_input_is_somatic
  par_sites_only_vcf_output
  par_disable_sequence_dictionary_validation
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

if [[ -n "$par_dbsnp_index" && -z "$par_dbsnp" ]]; then
  echo "Error: --dbsnp_index requires --dbsnp." >&2
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
tmp_dir=$(mktemp -d "$meta_temp_dir/gatk4_combinegvcfs.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

staged_reference=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict") || exit 1
# GATK reads a plain --variant VCF without an index, unless it must traverse
# by intervals (any interval argument)
variant_index_optional=$([[ "$(gatk_uses_intervals)" == "true" ]] && echo false || echo true)
stage_vcfs_with_indexes "$tmp_dir" "$par_variant" "$par_variant_index" "--variant" "--variant" variant_args "$variant_index_optional"
if [[ -n "$par_dbsnp" ]]; then
  # GATK reads a plain --dbsnp VCF without an index
  staged_dbsnp=$(stage_vcf_with_index "$tmp_dir" "$par_dbsnp" "$par_dbsnp_index" dbsnp true) || exit 1
fi
stage_interval_arguments "$tmp_dir" interval_args

# Convert semicolon-separated multi-value arguments into repeated flags
split_multiple_to_flags "$par_annotation" "--annotation" annotation_args
split_multiple_to_flags "$par_annotation_group" "--annotation-group" annotation_group_args
split_multiple_to_flags "$par_annotations_to_exclude" "--annotations-to-exclude" annotations_to_exclude_args
split_multiple_to_flags "$par_founder_id" "--founder-id" founder_id_args

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Build command arguments array
cmd_args=(
  ${par_allow_old_rms_mapping_quality_annotation_data:+--allow-old-rms-mapping-quality-annotation-data}
  "${annotation_args[@]}"
  "${annotation_group_args[@]}"
  "${annotations_to_exclude_args[@]}"
  ${par_break_bands_at_multiples_of:+--break-bands-at-multiples-of "$par_break_bands_at_multiples_of"}
  ${par_call_genotypes:+--call-genotypes}
  ${par_combine_variants_distance:+--combine-variants-distance "$par_combine_variants_distance"}
  ${par_convert_to_base_pair_resolution:+--convert-to-base-pair-resolution}
  ${staged_dbsnp:+--dbsnp "$staged_dbsnp"}
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
  ${par_sites_only_vcf_output:+--sites-only-vcf-output}
  ${par_create_output_variant_index:+--create-output-variant-index "$par_create_output_variant_index"}
  ${par_disable_sequence_dictionary_validation:+--disable-sequence-dictionary-validation}
  ${par_lenient:+--lenient}
)

# Run GATK CombineGVCFs
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" CombineGVCFs \
  "${variant_args[@]}" \
  --reference "$staged_reference" \
  --output "$par_output" \
  "${cmd_args[@]}" \
  --tmp-dir "$tmp_dir"

# Move the output index to --output_index, if requested
if [[ -n "$move_output_index" && ! "$gatk_output_index" -ef "$par_output_index" ]]; then
  mv "$gatk_output_index" "$par_output_index"
fi
