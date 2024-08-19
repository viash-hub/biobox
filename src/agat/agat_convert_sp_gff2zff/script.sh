#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

agat_convert_sp_gff2zff.pl \
  --gff "$par_gff" \
  --fasta "$par_fasta" \
  --output "$par_output" \
  ${par_config:+--config "${par_config}"} 

# move the output to the given directory
if [ -n "$par_output_dir" ]; then
  mkdir -p "$par_output_dir"
  mv "$par_output".ann "$par_output_dir"/"$par_output".ann
  mv "$par_output".dna "$par_output_dir"/"$par_output".dna 
fi