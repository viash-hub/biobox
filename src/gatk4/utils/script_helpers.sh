#!/bin/bash

# GATK4-specific runtime script helper functions for biobox components
#
# Source this file (alongside the component's own script.sh) via:
#   source "$meta_resources_dir/gatk4/script_helpers.sh"

# Split a semicolon-separated `multiple: true` argument value into a bash
# array of repeated "--flag value" pairs, e.g. "a;b;c" with flag
# "--read-filter" becomes ("--read-filter" "a" "--read-filter" "b"
# "--read-filter" "c"). Produces an empty array if the value is unset/empty.
#
# Usage: split_multiple_to_flags "$par_value" "--flag-name" result_array_name
split_multiple_to_flags() {
  local value="$1"
  local flag="$2"
  local -n result_ref="$3"

  result_ref=()
  if [[ -n "$value" ]]; then
    local items=()
    IFS=';' read -ra items <<< "$value"
    local item
    for item in "${items[@]}"; do
      result_ref+=("$flag" "$item")
    done
  fi
}

# Symlink a reference FASTA and its .fai/.dict into a temp dir under matching
# basenames, so GATK can find them
#
# Usage: staged=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
stage_reference_trio() {
  local tmp_dir="$1"
  local reference="$2"
  local reference_fai="$3"
  local reference_dict="$4"

  ln -s "$(readlink -f "$reference")" "$tmp_dir/reference.fasta"
  ln -s "$(readlink -f "$reference_fai")" "$tmp_dir/reference.fasta.fai"
  ln -s "$(readlink -f "$reference_dict")" "$tmp_dir/reference.dict"

  echo "$tmp_dir/reference.fasta"
}

# Symlink a BAM and its .bai companion into a temp dir under matching
# basenames, so GATK can find them
#
# Usage: staged=$(stage_bam_bai "$tmp_dir" "$par_input" "$par_bai")
stage_bam_bai() {
  local tmp_dir="$1"
  local input_bam="$2"
  local input_bai="$3"

  ln -s "$(readlink -f "$input_bam")" "$tmp_dir/sample.bam"
  ln -s "$(readlink -f "$input_bai")" "$tmp_dir/sample.bai"

  echo "$tmp_dir/sample.bam"
}
