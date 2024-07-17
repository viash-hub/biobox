#!/bin/bash

## VIASH START
## VIASH END

agat_convert_sp_gff2gtf.pl \
  -i "$par_gff" \
  -o "$par_output" \
  ${par_gtf_version:+--gtf_version "${par_gtf_version}"} \
  ${par_config:+--config "${par_config}"}
