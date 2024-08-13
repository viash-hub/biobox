#!/bin/bash

set -e

## VIASH START
## VIASH END

# create temporary directory
tmp_dir=$(mktemp -d -p "$meta_temp_dir" "${meta_functionality_name}_XXXXXX")
mkdir -p "$tmp_dir/temp"

# create detailed_results directory if variable is set and directory does not exist
if [[ ! -z "$par_detailed_results" ]] && [[ ! -d "$par_detailed_results" ]]; then
  mkdir -p "$par_detailed_results"
fi

# replace comma with semicolon
par_feature_type=$(echo $par_feature_type | tr ',' ';')
par_extra_attributes=$(echo $par_extra_attributes | tr ',' ';')

# unset flag variables
unset_if_false=(
    par_feature_level
    par_overlapping
    par_largest_overlap
    par_multi_mapping
    par_fraction
    par_split_only
    par_non_split_only
    par_primary
    par_ignore_dup
    par_paired
    par_count_read_pairs
    par_both_aligned
    par_check_pe_dist
    par_same_strand
    par_donotsort
    par_by_read_group
    par_long_reads
    par_verbose
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

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
  ${par_detailed_results:+--Rpath "${par_detailed_results}"} \
  ${par_detailed_results_format:+-R "${par_detailed_results_format}"} \
  ${par_max_M_op:+--maxMOp "${par_max_M_op}"} \
  ${par_verbose:+--verbose} \
  ${meta_cpus:+-T "${meta_cpus}"} \
  --tmpDir "$tmp_dir/temp" \
  -a "$par_annotation" \
  -o "$tmp_dir/output.txt" \
  "${input[*]}"

[[ ! -z "$par_counts" ]] && mv "$tmp_dir/output.txt" "$par_counts"
[[ ! -z "$par_summary" ]] && mv "$tmp_dir/output.txt.summary" "$par_summary"
if [[ ! -z "$par_junctions" ]] && [[ -e "$tmp_dir/output.txt.jcounts" ]]; then 
  mv "$tmp_dir/output.txt.jcounts" "$par_junctions"
fi
