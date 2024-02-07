#!/bin/bash

if [ -z $par_output ]; then
	par_output=.
else
	mkdir -p "$par_output"
fi

echo "par_pair_adapters: $par_pair_adapters"
echo "par_pair_filter: $par_pair_filter"
echo "par_interleaved: $par_interleaved"
echo "par_error_rate: $par_error_rate"
echo "par_no_indels: $par_no_indels"
echo "par_times: $par_times"
echo "par_overlap: $par_overlap"
echo "par_match_read_wildcards: $par_match_read_wildcards"
echo "no_match_adapter_wildcards: $no_match_adapter_wildcards"
echo "par_action: $par_action"
echo "par_revcomp: $par_revcomp"
echo "par_cut: $par_cut"
echo "par_cutR2: $par_cutR2"
echo "par_nextseq_trim: $par_nextseq_trim"
echo "par_quality_cutoff: $par_quality_cutoff"
echo "par_quality_cutoffR2: $par_quality_cutoffR2"
echo "par_quality_base: $par_quality_base"
echo "par_poly_a: $par_poly_a"
echo "par_length: $par_length"
echo "par_trim_n: $par_trim_n"
echo "par_length_tag: $par_length_tag"
echo "par_strip_suffix: $par_strip_suffix"
echo "par_prefix: $par_prefix"
echo "par_suffix: $par_suffix"
echo "par_rename: $par_rename"
echo "par_zero_cap: $par_zero_cap"
echo "par_minimum_length: $par_minimum_length"
echo "par_maximum_length: $par_maximum_length"
echo "par_max_n: $par_max_n"
echo "par_max_expected_errors: $par_max_expected_errors"
echo "par_max_average_error_rate: $par_max_average_error_rate"
echo "par_discard_trimmed: $par_discard_trimmed"
echo "par_discard_untrimmed: $par_discard_untrimmed"
echo "par_discard_casava: $par_discard_casava"
echo "par_report: $par_report"
echo "par_json: $par_json"
echo "par_output: $par_output"
echo "par_fasta: $par_fasta"
echo "par_info_file: $par_info_file"

# Do we get explicit adapter sequences or a FASTA file?
# Let the underlying tool deal with inconsistant states.
adapter_mode_R1=""
if [ ! -z "${par_adapter_fasta+set}" ]; then
  adapter_mode_R1="fasta"
else
  adapter_mode_R1="plain"
fi

front_mode_R1=""
if [ ! -z "${par_front_fasta+set}" ]; then
  front_mode_R1="fasta"
else
  front_mode_R1="plain"
fi

anywhere_mode_R1=""
if [ ! -z "${par_anywhere_fasta+set}" ]; then
  anywhere_mode_R1="fasta"
else
  anywhere_mode_R1="plain"
fi

adapter_mode_R2=""
if [ ! -z "${par_adapter_fastaR2+set}" ]; then
  adapter_mode_R2="fasta"
else
  adapter_mode_R2="plain"
fi

front_mode_R2=""
if [ ! -z "${par_front_fastaR2+set}" ]; then
  front_mode_R2="fasta"
else
  front_mode_R2="plain"
fi

anywhere_mode_R2=""
if [ ! -z "${par_anywhere_fastaR2+set}" ]; then
  anywhere_mode_R2="fasta"
else
  anywhere_mode_R2="plain"
fi

echo "Running cutadapt"
echo
echo "Adapter settings"
echo "----------------"
echo "Adapter Mode  R1 : $adapter_mode_R1"
echo "Front Mode    R1 : $front_mode_R1"
echo "Anywhere Mode R1 : $anywhere_mode_R1"
echo "Adapter Mode  R2 : $adapter_mode_R2"
echo "Front Mode    R2 : $front_mode_R2"
echo "Anywhere Mode R2 : $anywhere_mode_R2"
echo

# Adapter arguments
#   - paired and single-end
#   - string and fasta
###########################################################
echo ">> Parsing arguments dealing with adapters"
adapter_args=$(echo \
  ${par_adapter:+--adapter "${par_adapter}"} \
  ${par_adapter_fasta:+--adapter "file:${par_adapter_fasta}"} \
  ${par_front:+--front "${par_front}"} \
  ${par_front_fasta:+--front "file:${par_front_fasta}"} \
  ${par_anywhere:+--anywhere "${par_anywhere}"} \
  ${par_anywhere_fasta:+--anywhere "file:${par_anywhere_fasta}"} \
  ${par_adapterR2:+--adapterR2 "${par_adapterR2}"} \
  ${par_adapterR2_fasta:+--adapterR2 "file:${par_adapterR2_fasta}"} \
  ${par_frontR2:+--frontR2 "${par_frontR2}"} \
  ${par_frontR2_fasta:+--frontR2 "file:${par_frontR2_fasta}"} \
  ${par_anywhereR2:+--anywhereR2 "${par_anywhereR2}"} \
  ${par_anywhereR2_fasta:+--anywhereR2 "file:${par_anywhereR2_fasta}"}
)
echo "Arguments to cutadapt:"
echo "$adapter_args"
echo

# Paired-end options
###########################################################
echo ">> Parsing arguments for paired-end reads"
[[ "$par_pair_adapters" == "false" ]] && unset par_pair_adapters
[[ "$par_interleaved" == "false" ]] && unset par_interleaved

paired_args=$(echo \
  ${par_pair_adapters:+--pair-adapters} \
  ${par_pair_filter:+--pair-filter "${par_pair_filter}"} \
  ${par_interleaved:+--interleaved} \

)
echo "Arguments to cutadapt:"
echo $paired_args
echo

# Input arguments 
###########################################################
echo ">> Parsing input arguments"
[[ "$par_no_indels" == "true" ]] && unset par_no_indels
[[ "$par_match_read_wildcards" == "false" ]] && unset par_match_read_wildcards
[[ "$par_no_match_adapter_wildcards" == "true" ]] && unset par_no_match_adapter_wildcards
[[ "$par_revcomp" == "false" ]] && unset par_revcomp

input_args=$(echo \
  ${par_error_rate:+-error-rate "${par_error_rate}"} \
  ${par_no_indels:+--no-indels} \
  ${par_times:+--times "${par_times}"} \
  ${par_overlap:+--overlap "${par_overlap}"} \
  ${par_match_read_wildcards:+--match-read-wildcards} \
  ${par_no_match_adapter_wildcards:+--no-match-adapter-wildcards} \
  ${par_action:+--action "${par_action}"} \
  ${par_revcomp:+--revcomp} \
)
echo "Arguments to cutadapt:"
echo $input_args
echo

# Read modifications
###########################################################
echo ">> Parsing read modification arguments"
[[ "$par_poly_a" == "false" ]] && unset par_poly_a
[[ "$par_trim_n" == "false" ]] && unset par_trim_n
[[ "$par_zero_cap" == "false" ]] && unset par_zero_cap

mod_args=$(echo \
  ${par_cut:+--cut "${par_cut}"} \
  ${par_cutR2:+--cutR2 "${par_cutR2}"} \
  ${par_nextseq_trim:+--nextseq-trim "${par_nextseq_trim}"} \
  ${par_quality_cutoff:+--quality-cutoff "${par_quality_cutoff}"} \
  ${par_quality_cutoffR2:+--quality-cutoffR2 "${par_quality_cutoffR2}"} \
  ${par_quality_base:+--quality-base "${par_quality_base}"} \
  ${par_poly_a:+--poly-a} \
  ${par_length:+--length "${par_length}"} \
  ${par_trim_n:+--trim-n} \
  ${par_length_tag:+--length-tag "${par_length_tag}"} \
  ${par_strip_suffix:+--strip-suffix "${par_strip_suffix}"} \
  ${par_prefix:+--prefix "${par_prefix}"} \
  ${par_suffix:+--suffix "${par_suffix}"} \
  ${par_rename:+--rename "${par_rename}"} \
  ${par_zero_cap:+--zero-cap} \
)
echo "Arguments to cutadapt:"
echo $mod_args
echo

# Filtering of processed reads arguments
###########################################################
echo ">> Filtering of processed reads arguments"
[[ "$par_discard_trimmed" == "false" ]] && unset par_discard_trimmed
[[ "$par_discard_untrimmed" == "false" ]] && unset par_discard_untrimmed
[[ "$par_discard_casava" == "false" ]] && unset par_discard_casava

filter_args=$(echo \
  ${par_minimum_length:+--minimum-length "${par_minimum_length}"} \
  ${par_maximum_length:+--maximum-length "${par_maximum_length}"} \
  ${par_max_n:+--max-n "${par_max_n}"} \
  ${par_max_expected_errors:+--max-expected-errors "${par_max_expected_errors}"} \
  ${par_max_average_error_rate:+--max-average-error-rate "${par_max_average_error_rate}"} \
  ${par_discard_trimmed:+--discard-trimmed} \
  ${par_discard_untrimmed:+--discard-untrimmed} \
  ${par_discard_casava:+--discard-casava} \
)
echo "Arguments to cutadapt:"
echo $filter_args
echo

# Output arguments
# We write the output to a directory rather than
# individual files.
###########################################################
echo ">> Output arguments"
[[ "$par_json" == "false" ]] && unset par_json
[[ "$par_fasta" == "false" ]] && unset par_fasta
[[ "$par_info_file" == "false" ]] && unset par_info_file

# 	-o "$par_outputDir/{name}_R1_001.fastq" \
# 	-p "$par_outputDir/{name}_R2_001.fastq" \

output_args=$(echo \
  ${par_report:+--report "${par_report}"} \
  ${par_json:+--json} \
  -o "$par_output/{name}_R1_001.fastq" \
  -p "$par_output/{name}_R1_001.fastq" \
  ${par_fasta:+--fasta} \
  ${par_info_file:+--info-file} \
)
echo "Arguments to cutadapt:"
echo $output_args
echo

echo ">> Full CLI to be run:"
cli=$(echo "cutadapt" \
  $adapter_args \
  $paired_args \
  $input_args \
  $mod_args \
  $filter_args \
  $output_args
)

echo $cli

# $( "$cli" ) > $par_output/report.txt
