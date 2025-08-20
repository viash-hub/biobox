#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_counts" == "false" ]] && unset par_counts
[[ "$par_both" == "false" ]] && unset par_both
[[ "$par_strand" == "false" ]] && unset par_strand
[[ "$par_different_strand" == "false" ]] && unset par_different_strand

# Convert semicolon-separated files to array
IFS=';' read -ra files_array <<< "$par_files"

# Convert semicolon-separated names to array if provided
if [[ -n "${par_names}" ]]; then
    IFS=';' read -ra names_array <<< "$par_names"
fi

# Build command arguments array
cmd_args=(
    -i "$par_input"
    ${par_names:+-names "${names_array[@]}"}
    ${par_counts:+-counts}
    ${par_both:+-both}
    ${par_strand:+-s}
    ${par_different_strand:+-S}
    -files "${files_array[@]}"
)

# Execute bedtools annotate
bedtools annotate "${cmd_args[@]}" > "$par_output"
