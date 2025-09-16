#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters
[[ "$par_first_allele_only" == "false" ]] && unset par_first_allele_only
[[ "$par_split_by_id" == "false" ]] && unset par_split_by_id
[[ "$par_verbose" == "false" ]] && unset par_verbose

# Handle multiple input files (semicolon-separated)
IFS=';' read -ra input_files <<< "$par_input"

# Validate input files (maximum 2 allowed)
if [[ ${#input_files[@]} -gt 2 ]]; then
  echo "Error: Maximum of two input files allowed" >&2
  exit 1
fi

# Build command array
cmd_args=(
  "bcftools" "stats"
  ${par_af_bins:+--af-bins "$par_af_bins"}
  ${par_af_tag:+--af-tag "$par_af_tag"}
  ${par_collapse:+-c "$par_collapse"}
  ${par_depth:+-d "$par_depth"}
  ${par_exclude:+-e "$par_exclude"}
  ${par_exons:+-E "$par_exons"}
  ${par_apply_filters:+-f "$par_apply_filters"}
  ${par_fasta_ref:+-F "$par_fasta_ref"}
  ${par_first_allele_only:+--1st-allele-only}
  ${par_include:+-i "$par_include"}
  ${par_split_by_id:+-I}
  ${par_regions:+-r "$par_regions"}
  ${par_regions_file:+-R "$par_regions_file"}
  ${par_regions_overlap:+--regions-overlap "$par_regions_overlap"}
  ${par_samples:+-s "$par_samples"}
  ${par_samples_file:+-S "$par_samples_file"}
  ${par_targets:+-t "$par_targets"}
  ${par_targets_file:+-T "$par_targets_file"}
  ${par_targets_overlap:+--targets-overlap "$par_targets_overlap"}
  ${par_user_tstv:+-u "$par_user_tstv"}
  ${par_verbose:+-v}
)

# Add input files
for file in "${input_files[@]}"; do
  cmd_args+=("$file")
done

# Execute command and redirect output
"${cmd_args[@]}" > "$par_output"

