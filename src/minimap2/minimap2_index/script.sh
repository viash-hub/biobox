#!/usr/bin/env bash

set -euo pipefail

## VIASH START
## VIASH END

echo "Creating minimap2 index..."
echo "Input FASTA: $par_input"
echo "Output index: $par_output"

# Validate input parameters
if [ -z "${par_input:-}" ]; then
  echo "Error: Input file is required."
  usage
fi

if [ -z "${par_output:-}" ]; then
  echo "Error: Output file is required."
  usage
fi

# Check if input file exists
if [ ! -f "$par_input" ]; then
  echo "Error: Input file '$par_input' not found!"
  exit 1
fi

# Run minimap2 with the -d flag to create the index
minimap2 \
  -d "$par_output" \
  "$par_input"

# Check if the output index file was created successfully
if [ ! -f "$par_output" ]; then
  echo "Error: Index file '$par_output' was not created."
  exit 1
fi

echo "Minimap2 indexing complete."
