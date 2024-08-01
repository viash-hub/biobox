#!/bin/bash

## VIASH START
## VIASH END

agat_convert_sp_gxf2gxf.pl \
  -g "$par_gxf" \
  -o "$par_output" \
  ${par_config:+--config "${par_config}"}
