#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Build command arguments array
cmd_args=(
    -i "$par_input"
    -c "$par_columns"
)

# Execute bedtools expand
bedtools expand "${cmd_args[@]}" > "$par_output"
