#!/bin/bash

set -euo pipefail

## VIASH START
## VIASH END

# unset boolean flags that are false
[[ "${par_fasta:-false}" == "false" ]] && unset par_fasta

# building command arguments
cmd_args=(
  ${meta_cpus:+-p "$meta_cpus"}
  -x "$par_index_dir/$par_index_prefix"
  -S "$par_output_sam"
  ${par_fasta:+-f}
)

# single-end vs paired-end mode handling
if [[ -n "${par_input_r2:-}" ]]; then
  echo "**Running in paired-end mode**"
  cmd_args+=(-1 "$par_input" -2 "$par_input_r2")
else
  echo "**Running in single-end mode**"
  cmd_args+=(-U "$par_input")
fi

# run hisat2 align
hisat2 "${cmd_args[@]}"

echo "✓ Alignment complete"
echo "Please find HISAT2 alignment output here: $par_output_sam"