#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

agat_sp_add_splice_sites.pl \
  --gff "$par_gff" \
  --output "$par_output" \
  ${par_config:+--config "${par_config}"}
