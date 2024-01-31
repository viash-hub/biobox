#!/bin/bash

## VIASH START
## VIASH END

tmp_dir=$(mktemp -d -p "$meta_temp_dir" "${meta_functionality_name}_XXXXXXXXX")

if [[ $par_r_path ]] && [[ ! -d "$par_r_path" ]]; then
    mkdir -p "$par_r_path"
fi

[[ "$par_feature_level" == "false" ]] && unset par_feature_level
[[ "$par_overlapping" == "false" ]] && unset par_overlapping
[[ "$par_largest_overlap" == "false" ]] && unset par_largest_overlap
[[ "$par_multi_mapping" == "false" ]] && unset par_multi_mapping
[[ "$par_fraction" == "false" ]] && unset par_fraction
[[ "$par_split_only" == "false" ]] && unset par_split_only
[[ "$par_non_split_only" == "false" ]] && unset par_non_split_only
[[ "$par_primary" == "false" ]] && unset par_primary
[[ "$par_ignore_dup" == "false" ]] && unset par_ignore_dup
[[ "$par_junctions" == "false" ]] && unset par_junctions
[[ "$par_paired" == "false" ]] && unset par_paired
[[ "$par_count_read_pairs" == "false" ]] && unset par_count_read_pairs
[[ "$par_both_aligned" == "false" ]] && unset par_both_aligned
[[ "$par_check_pe_dist" == "false" ]] && unset par_check_pe_dist
[[ "$par_same_strand" == "false" ]] && unset par_same_strand
[[ "$par_donotsort" == "false" ]] && unset par_donotsort
[[ "$par_by_read_group" == "false" ]] && unset par_by_read_group
[[ "$par_long_reads" == "false" ]] && unset par_long_reads
[[ "$par_verbose" == "false" ]] && unset par_verbose
[[ "$par_version" == "false" ]] && unset par_version

IFS=";" read -ra input <<< $par_input

featureCounts \
  ${par_format:+-F "${par_format}"} \
  ${par_feature_type:+-t "${par_feature_type}"} \
  ${par_attribute_type:+-g "${par_attribute_type}"} \
  ${par_extra_attributes:+--extraAttributes "${extra_attributes}"} \
  ${par_chrom_alias:+-A "${par_chrom_alias}"} \
  ${par_feature_level:+-f} \
  ${par_overlapping:+-O} \
  ${par_min_overlap:+--minOverlap "${par_min_overlap}"} \
  ${par_frac_overlap:+--fracOverlap "${par_frac_overlap}"} \
  ${par_frac_overlap_feature:+--fracOverlapFeature "${par_frac_overlap_feature}"} \
  ${par_largest_overlap:+--largestOverlap} \
  ${par_non_overlap:+--nonOverlap "${par_non_overlap}"} \
  ${par_non_overlap_feature:+--nonOverlapFeature "${par_non_overlap_feature}"} \
  ${par_read_extension5:+--readExtension5 "${par_read_extension5}"} \
  ${par_read_extension3:+--readExtension3 "${par_read_extension3}"} \
  ${par_read2pos:+--read2pos "${par_read2pos}"} \
  ${par_multi_mapping:+-M} \
  ${par_fraction:+--fraction} \
  ${par_min_map_quality:+-Q "${par_min_map_quality}"} \
  ${par_split_only:+--splitOnly} \
  ${par_non_split_only:+--nonSplitOnly} \
  ${par_primary:+--primary} \
  ${par_ignore_dup:+--ignoreDup} \
  ${par_strand:+-s "${par_strand}"} \
  ${par_junctions:+-J} \
  ${par_ref_fasta:+-G "${par_ref_fasta}"} \
  ${par_paired:+-p} \
  ${par_count_read_pairs:+--countReadPairs} \
  ${par_both_aligned:+-B} \
  ${par_check_pe_dist:+-P} \
  ${par_min_length:+-d "${par_min_length}"} \
  ${par_max_length:+-D "${par_max_length}"} \
  ${par_same_strand:+-C} \
  ${par_donotsort:+--donotsort} \
  ${par_by_read_group:+--byReadGroup} \
  ${par_long_reads:+-L} \
  ${par_r_path:+--r_path "${par_r_path}"} \
  ${par_detailed_results_format:+-R "${par_detailed_results_format}"} \
  ${par_tmpdir:+--tmpDir "${par_tmpdir}"} \
  ${par_max_M_op:+--maxMOp "${par_max_M_op}"} \
  ${par_verbose:+--verbose} \
  ${meta_cpus:+-T "${meta_cpus}"} \
  -a "$par_annotation" \
  -o "$par_output" \
  "${input[*]}"

[[ ! -z "$par_output_summary" ]] && mv "$par_output.summary" "$par_output_summary"
[[ ! -z "$par_output_junctions" ]] && mv "$par_output.jcounts" "$par_output_junctions"