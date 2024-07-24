#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_casava" == "false" ]] && unset par_casava
[[ "$par_nano" == "false" ]] && unset par_nano
[[ "$par_nofilter" == "false" ]] && unset par_nofilter
[[ "$par_extract" == "false" ]] && unset par_extract
[[ "$par_noextract" == "false" ]] && unset par_noextract
[[ "$par_nogroup" == "false" ]] && unset par_nogroup
[[ "$par_quiet" == "false" ]] && unset par_quiet

# Create input array 
IFS=";" read -ra input <<< $par_input

run fastqc
fastqc \
  ${par_outdir:+--outdir "$par_outdir"} \
  ${par_casava:+--casava} \
  ${par_nano:+--nano} \
  ${par_nofilter:+--nofilter} \
  ${par_extract:+--extract} \
  ${par_java:+--java "$par_java"} \
  ${par_noextract:+--noextract} \
  ${par_nogroup:+--nogroup} \
  ${par_min_length:+--min_length "$par_min_length"} \
  ${par_format:+--format "$par_format"} \
  ${par_contaminants:+--contaminants "$par_contaminants"} \
  ${par_adapters:+--adapters "$par_adapters"} \
  ${par_limits:+--limits "$par_limits"} \
  ${par_kmers:+--kmers "$par_kmers"} \
  ${par_quiet:+--quiet} \
  ${meta_cpus:+--threads "$meta_cpus"} \
  ${meta_temp_dir:+--dir "$meta_temp_dir"} \
  ${par_input:+ ${input[*]}}

# Get input directory
input_dir=$(dirname ${input[0]})

# Both outputs args passed
if [[ -n "$par_html" ]] && [[ -n "$par_zip" ]]; then
  for i in "${!input[@]}"; do
    # Substitute .bam, .sam, .bam_mapped, .sam_mapped, and .fastq to .fq
    input_html=$(echo ${input[$i]} | sed 's/\(\.bam\|\.sam\|\.bam_mapped\|\.sam_mapped\|\.fastq\)/.fq/g') 
    sample_name=$(basename $input_html .fq)
    input_zip="$input_dir/${sample_name}_fastqc.zip"
    input_html="$input_dir/${sample_name}_fastqc.html"
    zip_file="${par_zip//\*/$sample_name}"
    html_file="${par_html//\*/$sample_name}"
    mv "$input_zip" "$zip_file"
    mv "$input_html" "$html_file"
  done
# Only html output arg passed
elif [[ -n "$par_html" ]]; then
  for i in "${!input[@]}"; do
    # Substitute .bam, .sam, .bam_mapped, .sam_mapped, and .fastq to .fq
    input_html=$(echo ${input[$i]} | sed 's/\(\.bam\|\.sam\|\.bam_mapped\|\.sam_mapped\|\.fastq\)/.fq/g') 
    sample_name=$(basename $input_html .fq)
    input_html="$input_dir/${sample_name}_fastqc.html"
    html_file="${par_html//\*/$sample_name}" # Substitute * with sample name
    mv "$input_html" "$html_file"
  done
  rm "$input_dir"/*.zip
# Only zip output arg passed
elif [[ -n "$par_zip" ]]; then
  for i in "${!input[@]}"; do
    # Substitute .bam, .sam, .bam_mapped, .sam_mapped, and .fastq to .fq
    input_zip=$(echo ${input[$i]} | sed 's/\(\.bam\|\.sam\|\.bam_mapped\|\.sam_mapped\|\.fastq\)/.fq/g') 
    sample_name=$(basename $input_zip .fq)
    input_zip="$input_dir/${sample_name}_fastqc.zip"
    zip_file="${par_zip//\*/$sample_name}" # Substitute * with sample name
    mv "$input_zip" "$zip_file"
  done
  rm "$input_dir"/*.html
fi

