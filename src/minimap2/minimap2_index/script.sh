#!/usr/bin/env bash

set -euo pipefail

## VIASH START
## VIASH END

echo "Creating minimap2 index..."
echo "Input FASTA: $par_input"
echo "Output index: $par_output"

minimap2 \
  -d "$par_output" \
  "$par_input"

echo "Minimap2 indexing complete."
