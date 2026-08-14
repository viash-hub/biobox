#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

[[ "$par_zero_based" == "false" ]] && unset par_zero_based
[[ "$par_force" == "false" ]] && unset par_force
[[ "$par_csi" == "false" ]] && unset par_csi

# tabix always writes the index next to its input, so stage the input into a
# writable work dir before indexing it
work_dir=$(mktemp -d "$meta_temp_dir/tabix_index.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT

input_name=$(basename "$par_input")
staged_input="$work_dir/$input_name"
ln -s "$(realpath "$par_input")" "$staged_input"

cmd_args=(
  tabix
  ${par_preset:+--preset "$par_preset"}
  ${par_sequence:+--sequence "$par_sequence"}
  ${par_begin:+--begin "$par_begin"}
  ${par_end:+--end "$par_end"}
  ${par_skip_lines:+--skip-lines "$par_skip_lines"}
  ${par_comment:+--comment "$par_comment"}
  ${par_zero_based:+--zero-based}
  ${par_force:+--force}
  ${par_csi:+--csi}
  ${par_csi:+--min-shift "$par_min_shift"}
  ${meta_cpus:+--threads "$meta_cpus"}
  "$staged_input"
)

"${cmd_args[@]}"

# The resulting index name depends on input type and --csi (.tbi/.csi/.bai/.crai).
# Find the new file tabix created rather than hardcoding a suffix.
index_file=$(find "$work_dir" -maxdepth 1 -type f ! -name "$input_name")

cp "$index_file" "$par_output_index"
