#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset flags
[[ "$par_low_memory" == "false" ]] && unset par_low_memory
[[ "$par_meth" == "false" ]] && unset par_meth

# Create output directory
mkdir -p "$par_output"

if [ -n "$par_prefix" ]; then
    output_prefix="$par_output/$par_prefix"
else
    # Determine the index prefix for the output
    index_basename=$(basename "$par_input" .fasta)
    index_basename=$(basename "$index_basename" .fa)
    index_basename=$(basename "$index_basename" .fna)
    
    output_prefix="$par_output/$index_basename"
fi

# Build the command
cmd_args=(
    ${par_random_seed:+-s "$par_random_seed"}
    ${par_sa_sample_rate:+-u "$par_sa_sample_rate"}
    ${par_low_memory:+-l}
    ${par_block_size:+-b "$par_block_size"}
    ${meta_cpus:+-t "$meta_cpus"}
    ${par_meth:+--meth}
    "$par_input"
    "$output_prefix"
)

# Run minibwa index
minibwa index "${cmd_args[@]}"
