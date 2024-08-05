#!/bin/bash

## VIASH START
## VIASH END

agat_convert_sp_gff2zff.pl \
  --gff "$par_gff" \
  --fasta "$par_fasta" \
  -output "$par_output" \
  ${par_config:+--config "${par_config}"} 
