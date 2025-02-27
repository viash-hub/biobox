#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

agat_sp_clipN_seqExtremities_and_fixCoordinates.pl \
  --gff "$par_gff" \
  --fasta "$par_fasta" \
  --of "$par_output_fasta" \
  --og "$par_output_gff" \
  ${par_config:+--config "${par_config}"}