#!/usr/bin/env bash

set -euo pipefail

## VIASH START
par_index_dir="/Users/cydricgeyskens/Documents/viashhub/hisat2/index"
par_index_prefix=genome
par_input="/Users/cydricgeyskens/Documents/viashhub/hisat2/reads/reads_1.fa"
par_input_r2="/Users/cydricgeyskens/Documents/viashhub/hisat2/reads/reads_2.fa"
par_output_sam="/Users/cydricgeyskens/Documents/viashhub/hisat2/align_output.sam"
par_fasta="true"
## VIASH END

# unset boolean flags that are false
[[ "${par_fasta:-false}" == "false" ]] && unset par_fasta

echo "Aligning reads with HISAT2..."
echo "Input directory: $par_index_dir"
echo "Index prefix: $par_index_prefix"

# building command arguments
cmd_args=(
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
echo "Please find HISAT2 alignemnt output here: $par_output_sam"