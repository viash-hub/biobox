#!/bin/bash

## VIASH START
## VIASH END

# exit on error
set -eo pipefail

# unset flags
unset_if_false=(
  par_casava
  par_nano
  par_nofilter
  par_extract
  par_noextract
  par_nogroup
  par_quiet
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Create input array 
IFS=";" read -ra input <<< $par_input

# Run fastqc
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

# Moves extracted output directory to current directory
if [[ -n "$par_extract" ]]; then
  mv "$input_dir/"*_fastqc $(pwd)
fi

# Both outputs args passed (html and zip)
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
else 
  mv "$input_dir"/*.html "$input_dir"/*.zip $(pwd)
fi

