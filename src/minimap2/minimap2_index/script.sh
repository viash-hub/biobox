#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_homopolymer_compressed" == "false" ]] && unset par_homopolymer_compressed

echo "Creating minimap2 index..."
echo "Input FASTA: $par_input"
echo "Output index: $par_output"

# minimap2 recommends giving -x first, so that the preset does not override the
# -k / -w / -H values that follow it.
cmd_args=(
  ${par_preset:+-x "$par_preset"}
  ${par_kmer_size:+-k "$par_kmer_size"}
  ${par_window_size:+-w "$par_window_size"}
  ${par_homopolymer_compressed:+-H}
  ${meta_cpus:+-t "$meta_cpus"}
  -d "$par_output"
  "$par_input"
)

minimap2 "${cmd_args[@]}"

echo "Minimap2 indexing complete."
