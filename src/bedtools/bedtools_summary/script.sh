#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Build command arguments array
cmd_args=(
    -i "$par_input"
    -g "$par_genome"
)

# Execute bedtools summary and redirect output to the specified output file
bedtools summary "${cmd_args[@]}" > "$par_output"
