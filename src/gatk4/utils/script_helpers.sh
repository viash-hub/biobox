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
