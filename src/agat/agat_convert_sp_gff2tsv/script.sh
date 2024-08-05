#!/bin/bash

## VIASH START
## VIASH END

agat_convert_sp_gff2tsv.pl \
  -gff "$par_gff" \
  -output "$par_output" \
  ${par_config:+--config "${par_config}"}
