#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset boolean flags that are false
[[ "$par_header" == "false" ]] && unset par_header

# Convert semicolon-separated files to array
IFS=';' read -ra files_array <<< "$par_files"

# Build command arguments array
cmd_args=(
  ${par_header:+-header}
  -i "${files_array[@]}"
)

# Execute bedtools unionbedg and redirect output to the specified output file
bedtools unionbedg "${cmd_args[@]}" > "$par_output"
