#!/usr/bin/env bash

set -euo pipefail

## VIASH START
par_input="/Users/cydricgeyskens/Documents/viashhub/hisat2/22_20-21M.fa"
par_index_dir="/Users/cydricgeyskens/Documents/viashhub/hisat2/index"
par_index_prefix=genome
## VIASH END

echo "Building hisat2 index ..."
echo "Input genome FASTA: $par_input"
echo "Output directory: $par_index_dir"
echo "Index prefix: $par_index_prefix"

mkdir -p "$par_index_dir"

# build command arguments
cmd_args=(
  "$par_input"
  "$par_index_dir/$par_index_prefix"
)

# run hisat2 build
hisat2-build "${cmd_args[@]}"

echo "✓ Building hisat2 index complete"