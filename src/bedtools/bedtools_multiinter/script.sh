#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_cluster" == "false" ]] && unset par_cluster
[[ "$par_header" == "false" ]] && unset par_header
[[ "$par_empty" == "false" ]] && unset par_empty

# Build command arguments array
cmd_args=(
    ${par_cluster:+--cluster}
    ${par_header:+--header}
    ${par_empty:+--empty}
    ${par_genome:+--genome "$par_genome"}
    ${par_filler:+--filler "$par_filler"}
)

# Handle multiple input files - Viash passes them as semicolon-separated string
IFS=';' read -ra input_files <<< "$par_input"
for file in "${input_files[@]}"; do
    cmd_args+=(-i "$file")
done

# Add names if provided
if [[ ${par_names+x} ]]; then
    IFS=';' read -ra names_array <<< "$par_names"
    cmd_args+=(--names "${names_array[@]}")
fi

# Execute bedtools multiinter
bedtools multiinter "${cmd_args[@]}" > "$par_output"
