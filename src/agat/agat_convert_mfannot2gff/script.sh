#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

agat_convert_mfannot2gff.pl \
  --mfannot "$par_mfannot" \
  --gff "$par_gff" \
  ${par_config:+--config "${par_config}"}
