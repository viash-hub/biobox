#!/bin/bash

set -eo pipefail

[[ ! -d $output_dir ]] && mkdir -p $par_output_dir

IFS=";" read -ra input <<< $par_input

unset_if_false=( par_phred33 par_phred64 par_fastqc par_illumina par_stranded_illumina par_nextera par_small_rna par_gzip par_dont_gzip par_no_report_file par_suppress_warn par_clock par_polyA par_rrbs par_non_directional par_keep par_paired par_retain_unpaired )

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

trim_galore \
    ${par_quality:+-q "${par_quality}"} \
    ${par_phred33:+--phred33} \
    ${par_phred64:+--phred64 } \
    ${par_fastqc:+--fastqc } \
    ${par_fastqc_args:+--fastqc_args "${par_fastqc_args}"} \
    ${par_adapter:+-a "${par_adapter}"} \
    ${par_adapter2:+-a2 "${par_adapter2}"} \
    ${par_illumina:+--illumina} \
    ${par_stranded_illumina:+--stranded_illumina} \
    ${par_nextera:+--nextera} \
    ${par_small_rna:+--small_rna} \
    ${par_consider_already_trimmed:+--consider_already_trimmed "${par_consider_already_trimmed}"} \
    ${par_max_length:+--max_length "${par_max_length}"} \
    ${par_stringency:+--stringency "${par_stringency}"} \
    ${par_error_rate:+-e "${par_error_rate}"} \
    ${par_gzip:+--gzip} \
    ${par_dont_gzip:+--dont_gzip} \
    ${par_length:+--length "${par_length}"} \
    ${par_max_n:+--max_n "${par_max_n}"} \
    ${par_trim_n:+--trim-n "${par_trim_n}"} \
    ${par_no_report_file:+--no_report_file} \
    ${par_suppress_warn:+--suppress_warn} \
    ${par_clip_R1:+--clip_R1 "${par_clip_R1}"} \
    ${par_clip_R2:+--clip_R2 "${par_clip_R2}"} \
    ${par_three_prime_clip_R1:+--three_prime_clip_R1 "${par_three_prime_clip_R1}"} \
    ${par_three_prime_clip_R2:+--three_prime_clip_R2 "${par_three_prime_clip_R2}"} \
    ${par_nextseq:+--nextseq "${par_nextseq}"} \
    ${par_basename:+-basename "${par_basename}"} \
    ${par_hardtrim5:+--hardtrim5 "${par_hardtrim5}"} \
    ${par_hardtrim3:+--hardtrim3 "${par_hardtrim3}"} \
    ${par_clock:+--clock} \
    ${par_polyA:+--polyA} \
    ${par_implicon:+--implicon "${par_implicon}"} \
    ${par_rrbs:+--rrbs} \
    ${par_non_directional:+--non_directional} \
    ${par_keep:+--keep} \
    ${par_paired:+--paired} \
    ${par_retain_unpaired:+--retain_unpaired} \
    ${par_length_1:+-r1 "${par_length_1}"} \
    ${par_length_2:+-r2 "${par_length_2}"} \
    ${par_cores:+-j "${par_cores}"} \
    -o $par_output_dir \
    ${input[*]}

# Split output directory to separate files
[[ ! -z "$par_trimmed_r1" ]] && mv "$par_output_dir/output.txt" "$trimmed_r1"
[[ ! -z "$par_trimmed_r2" ]] && mv "$par_output_dir/output.txt.summary" "$par_trimmed_r2"
[[ ! -z "$par_trimming_report_r1" ]] && mv "$par_output_dir/output.txt" "$par_trimming_report_r1"
[[ ! -z "$par_trimming_report_r2" ]] && mv "$par_output_dir/output.txt.summary" "$par_trimming_report_r2"
[[ ! -z "$par_trimmed_fastqc_html_1" ]] && mv "$par_output_dir/output.txt" "$par_trimmed_fastqc_html_1"
[[ ! -z "$par_trimmed_fastqc_html_2" ]] && mv "$par_output_dir/output.txt.summary" "$par_trimmed_fastqc_html_2"
[[ ! -z "$par_trimmed_fastqc_zip_1" ]] && mv "$par_output_dir/output.txt" "$par_trimmed_fastqc_zip_1"
[[ ! -z "$par_trimmed_fastqc_zip_2" ]] && mv "$par_output_dir/output.txt.summary" "$par_trimmed_fastqc_zip_2"

if [ $par_paired == "true" ]; then     
    input_r1=$(basename -- "${input[0]}")
    input_r2=$(basename -- "${input[1]}")
    [[ ! -z "$par_trimmed_r1" ]] && mv "$par_output_dir/*val_1*.f*q.*" "$par_trimmed_r1"
    [[ ! -z "$par_trimmed_r2" ]] && mv "$par_output_dir/*val_2*.f*q.*" "$par_trimmed_r2"
    [[ ! -z "$par_trimming_report_r1" ]] && mv "$par_output_dir/${input_r1}_trimming_report.txt" "$par_trimming_report_r1"
    [[ ! -z "$par_trimming_report_r2" ]] && mv "$par_output_dir/${input_r1}_trimming_report.txt" "$par_trimming_report_r2"
    [[ ! -z "$par_trimmed_fastqc_html_1" ]] && mv "$par_output_dir/*val_1*.html" "$par_trimmed_fastqc_html_1"
    [[ ! -z "$par_trimmed_fastqc_html_2" ]] && mv "$par_output_dir/*val_2*.html" "$par_trimmed_fastqc_html_2"
    [[ ! -z "$par_trimmed_fastqc_zip_1" ]] && mv "$par_output_dir/*val_1*.zip" "$par_trimmed_fastqc_zip_1"
    [[ ! -z "$par_trimmed_fastqc_zip_2" ]] && mv "$par_output_dir/*val_2*.zip" "$par_trimmed_fastqc_zip_2"
    [[ ! -z "$par_unpaired_r1" ]] && mv "$par_output_dir/*.unpaired_1.f*q.*" "$par_unpaired_r1"
    [[ ! -z "$par_unpaired_r2" ]] && mv "$par_output_dir/*.unpaired_2.f*q.*" "$par_unpaired_r2"
else
    input_r1=$(basename -- "${input[0]}")
    [[ ! -z "$par_trimmed_r1" ]] && mv "$par_output_dir/*_trimmed.fq*" "$trimmed_r1"
    [[ ! -z "$par_trimming_report_r1" ]] && mv "$par_output_dir/${input_r1}_trimming_report.txt" "$par_trimming_report_r1"
    [[ ! -z "$par_trimmed_fastqc_html_1" ]] && mv "$par_output_dir/*_trimmed_fastqc.html" "$par_trimmed_fastqc_html_1"
    [[ ! -z "$par_trimmed_fastqc_zip_1" ]] && mv "$par_output_dir/*_trimmed_fastqc.zip" "$par_trimmed_fastqc_zip_1"
