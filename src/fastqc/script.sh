#!/bin/bash

## VIASH START
## VIASH END

# exit on error
set -eo pipefail

# Check if both outputs are empty, at least one must be passed.
if [[ -z "$par_outdir" ]] && [[ -z "$par_html" ]] && [[ -z "$par_zip" ]] && [[ -z "$par_summary" ]] && [[ -z "$par_data" ]]; then
  echo "Error: At least one of the output arguments (--outdir, --html, --zip, --summary, and --data) must be passed."
  exit 1
fi

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

tmpdir=$(mktemp -d "${meta_temp_dir}/${meta_name}-XXXXXXXX")
function clean_up {
  rm -rf "$tmpdir"
}
trap clean_up EXIT

# Set output directory
if [[ -n "$par_outdir" ]]; then
  if [[ ! -d "$par_outdir" ]]; then
    mkdir -p "$par_outdir"
  fi
  output_dir="$par_outdir"
else
  output_dir="$tmpdir"
fi

# Create input array 
IFS=";" read -ra input <<< $par_input

# Run fastqc
fastqc \
  --extract \
  ${par_casava:+--casava} \
  ${par_nano:+--nano} \
  ${par_nofilter:+--nofilter} \
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
  --outdir "${output_dir}" \
  "${input[@]}"
 

# Move output files
for file in "${input[@]}"; do
  # Removes everything after the first dot of the basename
  sample_name=$(basename "${file}" | sed 's/\..*$//')
  if [[ -n "$par_html" ]]; then
    input_html="${output_dir}/${sample_name}_fastqc.html"
    if [[ ! -f "$input_html" ]]; then
      echo "WARNING: HTML file '$input_html' does not exist"
    else
      html_file="${par_html//\*/$sample_name}"
      cp "$input_html" "$html_file"
    fi
  fi
  if [[ -n "$par_zip" ]]; then
    input_zip="${output_dir}/${sample_name}_fastqc.zip"
    if [[ ! -f "$input_zip" ]]; then
      echo "WARNING: ZIP file '$input_zip' does not exist"
    else
      zip_file="${par_zip//\*/$sample_name}"
      cp "$input_zip" "$zip_file"
    fi
  fi
  if [[ -n "$par_summary" ]]; then
    summary_file="${output_dir}/${sample_name}_fastqc/summary.txt"
    if [[ ! -f "$summary_file" ]]; then
      echo "WARNING: Summary file '$summary_file' does not exist"
    else
      new_summary="${par_summary//\*/$sample_name}"
      cp "$summary_file" "$new_summary"
    fi
  fi
  if [[ -n "$par_data" ]]; then
    data_file="${output_dir}/${sample_name}_fastqc/fastqc_data.txt"
    if [[ ! -f "$data_file" ]]; then
      echo "WARNING: Data file '$data_file' does not exist"
    else
      new_data="${par_data//\*/$sample_name}"
      cp "$data_file" "$new_data"
    fi
  fi
done


