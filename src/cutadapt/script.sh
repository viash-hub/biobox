#!/bin/bash

## VIASH START
par_adapter='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;GGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
par_input='src/cutadapt/test_data/se/a.fastq'
par_report='full'
par_json='false'
par_output='output'
par_fasta='false'
par_info_file='false'
## VIASH END

if [ -z $par_output ]; then
	par_output=.
else
	mkdir -p "$par_output"
fi


# Init
###########################################################
echo "Running cutadapt"
echo
echo ">> Paired-end data or not?"

mode=""
if [[ -z $par_input_r2 ]]; then
  mode="se"
  echo "  Single end"
  input="$par_input"
else
  echo "  Paired end"
  mode="pe"
  input="$par_input $par_input_r2"
fi

# Adapter arguments
#   - paired and single-end
#   - string and fasta
###########################################################

function add_flags {
  local arg=$1
  local flag=$2
  local prefix=$3
  [[ -z $prefix ]] && prefix=""

  # This function should not be called if the input is empty
  # but check for it just in case
  if [[ -z $arg ]]; then
    return
  fi

  local output=""
  IFS=';' read -r -a array <<< "$arg"
  for a in "${array[@]}"; do
    output="$output $flag $prefix$a"
  done
  echo $output
}

echo ">> Parsing arguments dealing with adapters"
adapter_args=$(echo \
  ${par_adapter:+$(add_flags "$par_adapter" "--adapter")} \
  ${par_adapter_fasta:+$(add_flags "$par_adapter_fasta" "--adapter" "file:")} \
  ${par_front:+$(add_flags "$par_front" "--front")} \
  ${par_front_fasta:+$(add_flags "$par_front_fasta" "--front" "file:")} \
  ${par_anywhere:+$(add_flags "$par_anywhere" "--anywhere")} \
  ${par_anywhere_fasta:+$(add_flags "$par_anywhere_fasta" "--anywhere" "file:")} \

  ${par_adapter_r2:+$(add_flags "$par_adapter_r2" "--adapter_r2")} \
  ${par_adapter_fasta_r2:+$(add_flags "$par_adapter_fasta_r2" "--adapter_r2" "file:")} \
  ${par_front_r2:+$(add_flags "$par_front_r2" "--front_r2")} \
  ${par_front_fasta_r2:+$(add_flags "$par_front_fasta_r2" "--front_r2" "file:")} \
  ${par_anywhere_r2:+$(add_flags "$par_anywhere_r2" "--anywhere_r2")} \
  ${par_anywhere_fasta_r2:+$(add_flags "$par_anywhere_fasta_r2" "--anywhere_r2" "file:")} \
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
  ${par_interleaved:+--interleaved}
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
  ${par_error_rate:+--error-rate "${par_error_rate}"} \
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
  ${par_cut_r2:+--cut_r2 "${par_cut_r2}"} \
  ${par_nextseq_trim:+--nextseq-trim "${par_nextseq_trim}"} \
  ${par_quality_cutoff:+--quality-cutoff "${par_quality_cutoff}"} \
  ${par_quality_cutoff_r2:+--quality-cutoff_r2 "${par_quality_cutoff_r2}"} \
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

# Parse and transform the minimum and maximum length arguments
[[ -z $par_minimum_length   ]]

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

if [ $mode = "se" ]; then
  output_args=$(echo \
    ${par_report:+--report "${par_report}"} \
    ${par_json:+--json "${par_output}/report.json"} \
    --output "$par_output/{name}_R1_001.fastq" \
    ${par_fasta:+--fasta} \
    ${par_info_file:+--info-file} \
  )
else
  output_args=$(echo \
    ${par_report:+--report "${par_report}"} \
    ${par_json:+--json "${par_output}/report.json"} \
    --output "$par_output/{name}_R1_001.fastq" \
    --paired-output "$par_output/{name}_R2_001.fastq" \
    ${par_fasta:+--fasta} \
    ${par_info_file:+--info-file} \
  )
fi
echo "Arguments to cutadapt:"
echo $output_args
echo

# Full CLI
# Set the --cores argument to 0 unless meta_cpus is set
###########################################################
echo ">> Full CLI to be run:"

par_cpus=0
[[ ! -z $meta_cpus ]] && par_cpus=$meta_cpus

cli=$(echo \
  $input \
  $adapter_args \
  $paired_args \
  $input_args \
  $mod_args \
  $filter_args \
  $output_args \
  --cores $par_cpus
)

echo cutadapt $cli | sed -e 's/--/\r\n  --/g'

cutadapt $cli | tee $par_output/report.txt
