#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

agat_sp_add_introns.pl \
  -f "$par_gff" \
  -o "$par_output" \
  ${par_config:+--config "${par_config}"}
