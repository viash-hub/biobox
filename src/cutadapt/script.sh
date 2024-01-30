#!/bin/bash

if [ -z $par_output ]; then
	par_output=.
else
	mkdir -p "$par_output"
fi

echo "par_adapter: $par_adapter"
echo "par_front: $par_front"
echo "par_anywhere: $par_anywhere"
echo "par_adapter_fasta: $par_adapter_fasta"
echo "par_front_fasta: $par_front_fasta"
echo "par_anywhere_fasta: $par_anywhere_fasta"
echo "par_adapterR2: $par_adapterR2"
echo "par_frontR2: $par_frontR2"
echo "par_anywhereR2: $par_anywhereR2"
echo "par_adapterR2_fasta: $par_adapterR2_fasta"
echo "par_frontR2_fasta: $par_frontR2_fasta"
echo "par_anywhereR2_fasta: $par_anywhereR2_fasta"
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

# for f in $par_input; do
# 	[ ! -f "$f" ] && echo "The input file $f does not exist" && exit 1
# done
#
# barcodesFasta="barcodes.fasta"
#
# awk '{print ">"$1"\n""^"$1}' $par_barcodesFile >$barcodesFasta
#
# fastqFiles=$(echo $par_input | tr " " "\n")
# for file in $fastqFiles; do
# 	if echo "$file" | grep -q R1; then
# 		input_R1=$(echo $file | grep R1)
# 	fi
# 	if echo "$file" | grep -q R2; then
# 		input_R2=$(echo $file | grep R2)
# 	fi
# done
# demuxFilesIn="$input_R1 $input_R2"
#
# # Note to self:
# #   The eval is here to expand shell globs, this way it is possible to use
# #   for instance pointers to ".../...R?....fastq", but please use the double
# #   quotes and an absolute path!
# eval /usr/local/bin/cutadapt \
# 	-e "$par_e" \
# 	--no-indels \
# 	--action=none \
# 	--cores=0 \
# 	-g "file:$barcodesFasta" \
# 	-o "$par_outputDir/{name}_R1_001.fastq" \
# 	-p "$par_outputDir/{name}_R2_001.fastq" \
# 	"$demuxFilesIn" >"$par_report"
