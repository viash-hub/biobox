#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_use_64bit_names" == "false" ]] && unset par_use_64bit_names

# Create output directory if it doesn't exist
output_dir=$(dirname "$par_output_amb")
mkdir -p "$output_dir"

# Build the command
cmd_args=(
    ${par_algorithm:+-a "$par_algorithm"}
    ${par_prefix:+-p "$par_prefix"}
    ${par_block_size:+-b "$par_block_size"}
    ${par_use_64bit_names:+-6}
    "$par_input"
)

# Run bwa index
bwa index "${cmd_args[@]}"

# Determine the index base name based on what BWA actually created
if [ -n "$par_prefix" ]; then
    index_base="$par_prefix"
else
    index_base="$par_input"
fi

# Handle 64-bit naming convention
if [ "$par_use_64bit_names" = "true" ]; then
    suffix=".64"
else
    suffix=""
fi

# Copy index files to specified output locations (only if different)
if [ "${index_base}${suffix}.amb" != "$par_output_amb" ]; then
    cp "${index_base}${suffix}.amb" "$par_output_amb"
fi

if [ "${index_base}${suffix}.ann" != "$par_output_ann" ]; then
    cp "${index_base}${suffix}.ann" "$par_output_ann" 
fi

if [ "${index_base}${suffix}.bwt" != "$par_output_bwt" ]; then
    cp "${index_base}${suffix}.bwt" "$par_output_bwt"
fi

if [ "${index_base}${suffix}.pac" != "$par_output_pac" ]; then
    cp "${index_base}${suffix}.pac" "$par_output_pac"
fi

if [ "${index_base}${suffix}.sa" != "$par_output_sa" ]; then
    cp "${index_base}${suffix}.sa" "$par_output_sa"
fi
