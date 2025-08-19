#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_use_64bit_names" == "false" ]] && unset par_use_64bit_names

# Create output directory
mkdir -p "$par_output"

# Determine the index base name for the output
index_basename=$(basename "$par_input" .fasta)
index_basename=$(basename "$index_basename" .fa)
index_basename=$(basename "$index_basename" .fna)

# Set prefix to write directly to output directory
if [ -n "$par_prefix" ]; then
    # Use custom prefix in output directory
    output_prefix="$par_output/$par_prefix"
else
    # Use input filename (without extension) as prefix in output directory
    output_prefix="$par_output/$index_basename"
fi

# Build the command
cmd_args=(
    ${par_algorithm:+-a "$par_algorithm"}
    -p "$output_prefix"
    ${par_block_size:+-b "$par_block_size"}
    ${par_use_64bit_names:+-6}
    "$par_input"
)

# Run bwa index
bwa index "${cmd_args[@]}"
