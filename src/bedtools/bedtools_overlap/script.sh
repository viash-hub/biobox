#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Build command arguments array
cmd_args=(
  -i "$par_input"
  -cols "$par_cols"
)

# Execute bedtools overlap
bedtools overlap "${cmd_args[@]}" > "$par_output"
