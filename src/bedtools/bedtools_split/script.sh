#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Create output directory if it doesn't exist
mkdir -p "$par_output_dir"

# Build command arguments array
cmd_args=(
  -i "$par_input"
  -n "$par_number"
  ${par_prefix:+-p "$par_prefix"}
  ${par_algorithm:+-a "$par_algorithm"}
)

# Change to output directory and execute bedtools split
cd "$par_output_dir"
bedtools split "${cmd_args[@]}"
