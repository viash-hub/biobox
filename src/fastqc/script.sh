#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_casava" == "false" ]] && unset par_casava
[[ "$par_nano" == "false" ]] && unset par_nano
[[ "$par_nofilter" == "false" ]] && unset par_nofilter
[[ "$par_extract" == "false" ]] && unset par_extract
#[[ -z "$par_java" ]] && unset par_java
[[ "$par_noextract" == "false" ]] && unset par_noextract
[[ "$par_nogroup" == "false" ]] && unset par_nogroup
#[[ "$par_min_length" == "false" ]] && unset par_min_length
#[[ "$par_format" == "false" ]] && unset par_format
#[[ "$par_threads" == "false" ]] && unset par_threads
#[[ "$par_contaminants" == "false" ]] && unset par_contaminants
#[[ "$par_adapters" == "false" ]] && unset par_adapters
#[[ "$par_limits" == "false" ]] && unset par_limits
#[[ "$par_kmers" == "false" ]] && unset par_kmers
[[ "$par_quiet" == "false" ]] && unset par_quiet
#[[ "$par_dir" == "false" ]] && unset par_dir

# run fastqc
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
  ${par_threads:+--threads "$par_threads"} \
  ${par_contaminants:+--contaminants "$par_contaminants"} \
  ${par_adapters:+--adapters "$par_adapters"} \
  ${par_limits:+--limits "$par_limits"} \
  ${par_kmers:+--kmers "$par_kmers"} \
  ${par_quiet:+--quiet} \
  ${par_dir:+--dir "$par_dir"} \
  $par_input


# Retrives the directory of the input file
input_dir=$(dirname "$par_input")

#echo "input is: $par_input"

# Both outputs args passed
if [[ -n "$par_html" ]] && [[ -n "$par_zip" ]]; then
  mv "$input_dir"/*.html "$par_html"
  mv "$input_dir"/*.zip "$par_zip"
# Only html output arg passed
elif [[ -n "$par_html" ]]; then
  mv "$input_dir"/*.html "$par_html"
  rm "$input_dir"/*.zip
# Only zip output arg passed
elif [[ -n "$par_zip" ]]; then
  mv "$input_dir"/*.zip "$par_zip"
  rm "$input_dir"/*.html
fi


# Questions:
# Should I unzip the zip file and make multiple other outputs for the viash component?
# TODO: handle the output args if multiple files are passed



  