#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

agat_convert_sp_gff2bed.pl \
  --gff "$par_gff" \
  --output "$par_output" \
  ${par_nc:+--nc "${par_nc}"} \
  ${par_sub:+--sub "${par_sub}"} \
  ${par_config:+--config "${par_config}"}
