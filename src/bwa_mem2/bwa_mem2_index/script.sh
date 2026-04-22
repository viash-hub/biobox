#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

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
    -p "$output_prefix"
    "$par_input"
)

# Run bwa index
bwa-mem2 index "${cmd_args[@]}"
