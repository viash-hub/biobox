#!/bin/bash

set -euo pipefail

## VIASH START
## VIASH END

echo "Building hisat2 index ..."
echo "Input genome FASTA: $par_reference"
echo "Output directory: $par_index_dir"
echo "Index prefix: $par_index_prefix"

mkdir -p "$par_index_dir"

# build command arguments
cmd_args=(
  ${meta_cpus:+-p "$meta_cpus"}
  "$par_reference"
  "$par_index_dir/$par_index_prefix"
)

# run hisat2 build
hisat2-build "${cmd_args[@]}"

echo "✓ Building hisat2 index complete"