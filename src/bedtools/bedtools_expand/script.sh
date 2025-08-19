#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Add columns parameter (required)
args=(
  -i "$par_input"
  -c "$par_columns"
)

# Execute bedtools expand
bedtools expand "${args[@]}" > "$par_output"
