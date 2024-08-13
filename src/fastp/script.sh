#!/bin/bash

## VIASH START
## VIASH END

# disable flags
unset_if_false=(
    par_disable_adapter_trimming
    par_detect_adapter_for_pe
    par_merge
    par_include_unmerged
    par_interleaved_in
    par_fix_mgi_id
    par_phred64
    par_dont_overwrite
    par_verbose
    par_dedup
    par_dont_eval_duplication
    par_trim_poly_g
    par_disable_trim_poly_g
    par_trim_poly_x
    par_disable_quality_filtering
    par_disable_length_filtering
    par_low_complexity_filter
    par_umi
    par_overrepresentation_analysis
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# run command
fastp \
  -i "$par_in1" \
  -o "$par_out1" \
  ${par_in2:+--in2 "${par_in2}"} \
  ${par_out2:+--out2 "${par_out2}"} \
  ${par_unpaired1:+--unpaired1 "${par_unpaired1}"} \
  ${par_unpaired2:+--unpaired2 "${par_unpaired2}"} \
  ${par_failed_out:+--failed_out "${par_failed_out}"} \
  ${par_overlapped_out:+--overlapped_out "${par_overlapped_out}"} \
  ${par_json:+--json "${par_json}"} \
  ${par_html:+--html "${par_html}"} \
  ${par_report_title:+--report_title "${par_report_title}"} \
  ${par_disable_adapter_trimming:+--disable_adapter_trimming} \
  ${par_detect_adapter_for_pe:+--detect_adapter_for_pe} \
  ${par_adapter_sequence:+--adapter_sequence "${par_adapter_sequence}"} \
  ${par_adapter_sequence_r2:+--adapter_sequence_r2 "${par_adapter_sequence_r2}"} \
  ${par_adapter_fasta:+--adapter_fasta "${par_adapter_fasta}"} \
  ${par_trim_front1:+--trim_front1 "${par_trim_front1}"} \
  ${par_trim_tail1:+--trim_tail1 "${par_trim_tail1}"} \
  ${par_max_len1:+--max_len1 "${par_max_len1}"} \
  ${par_trim_front2:+--trim_front2 "${par_trim_front2}"} \
  ${par_trim_tail2:+--trim_tail2 "${par_trim_tail2}"} \
  ${par_max_len2:+--max_len2 "${par_max_len2}"} \
  ${par_merge:+--merge} \
  ${par_merged_out:+--merged_out "${par_merged_out}"} \
  ${par_include_unmerged:+--include_unmerged} \
  ${par_interleaved_in:+--interleaved_in} \
  ${par_fix_mgi_id:+--fix_mgi_id} \
  ${par_phred64:+--phred64} \
  ${par_compression:+--compression "${par_compression}"} \
  ${par_dont_overwrite:+--dont_overwrite} \
  ${par_verbose:+--verbose} \
  ${par_reads_to_process:+--reads_to_process "${par_reads_to_process}"} \
  ${par_dedup:+--dedup} \
  ${par_dup_calc_accuracy:+--dup_calc_accuracy "${par_dup_calc_accuracy}"} \
  ${par_dont_eval_duplication:+--dont_eval_duplication} \
  ${par_trim_poly_g:+--trim_poly_g} \
  ${par_poly_g_min_len:+--poly_g_min_len "${par_poly_g_min_len}"} \
  ${par_disable_trim_poly_g:+--disable_trim_poly_g} \
  ${par_trim_poly_x:+--trim_poly_x} \
  ${par_poly_x_min_len:+--poly_x_min_len "${par_poly_x_min_len}"} \
  ${par_cut_front:+--cut_front "${par_cut_front}"} \
  ${par_cut_tail:+--cut_tail "${par_cut_tail}"} \
  ${par_cut_right:+--cut_right "${par_cut_right}"} \
  ${par_cut_window_size:+--cut_window_size "${par_cut_window_size}"} \
  ${par_cut_mean_quality:+--cut_mean_quality "${par_cut_mean_quality}"} \
  ${par_cut_front_window_size:+--cut_front_window_size "${par_cut_front_window_size}"} \
  ${par_cut_front_mean_quality:+--cut_front_mean_quality "${par_cut_front_mean_quality}"} \
  ${par_cut_tail_window_size:+--cut_tail_window_size "${par_cut_tail_window_size}"} \
  ${par_cut_tail_mean_quality:+--cut_tail_mean_quality "${par_cut_tail_mean_quality}"} \
  ${par_cut_right_window_size:+--cut_right_window_size "${par_cut_right_window_size}"} \
  ${par_cut_right_mean_quality:+--cut_right_mean_quality "${par_cut_right_mean_quality}"} \
  ${par_disable_quality_filtering:+--disable_quality_filtering} \
  ${par_qualified_quality_phred:+--qualified_quality_phred "${par_qualified_quality_phred}"} \
  ${par_unqualified_percent_limit:+--unqualified_percent_limit "${par_unqualified_percent_limit}"} \
  ${par_n_base_limit:+--n_base_limit "${par_n_base_limit}"} \
  ${par_average_qual:+--average_qual "${par_average_qual}"} \
  ${par_disable_length_filtering:+--disable_length_filtering} \
  ${par_length_required:+--length_required "${par_length_required}"} \
  ${par_length_limit:+--length_limit "${par_length_limit}"} \
  ${par_low_complexity_filter:+--low_complexity_filter} \
  ${par_complexity_threshold:+--complexity_threshold "${par_complexity_threshold}"} \
  ${par_filter_by_index1:+--filter_by_index1 "${par_filter_by_index1}"} \
  ${par_filter_by_index2:+--filter_by_index2 "${par_filter_by_index2}"} \
  ${par_filter_by_index_threshold:+--filter_by_index_threshold "${par_filter_by_index_threshold}"} \
  ${par_correction:+--correction} \
  ${par_overlap_len_require:+--overlap_len_require "${par_overlap_len_require}"} \
  ${par_overlap_diff_limit:+--overlap_diff_limit "${par_overlap_diff_limit}"} \
  ${par_overlap_diff_percent_limit:+--overlap_diff_percent_limit "${par_overlap_diff_percent_limit}"} \
  ${par_umi:+--umi} \
  ${par_umi_loc:+--umi_loc "${par_umi_loc}"} \
  ${par_umi_len:+--umi_len "${par_umi_len}"} \
  ${par_umi_prefix:+--umi_prefix "${par_umi_prefix}"} \
  ${par_umi_skip:+--umi_skip "${par_umi_skip}"} \
  ${par_umi_delim:+--umi_delim "${par_umi_delim}"} \
  ${par_overrepresentation_analysis:+--overrepresentation_analysis} \
  ${par_overrepresentation_sampling:+--overrepresentation_sampling "${par_overrepresentation_sampling}"} \
  ${meta_cpus:+--thread "${meta_cpus}"}
