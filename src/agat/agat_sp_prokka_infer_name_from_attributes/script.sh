#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# unset flags
[[ "$par_force" == "false" ]] && unset par_force

# run agat_sp_Prokka_inferNameFromAttributes.pl
agat_sp_Prokka_inferNameFromAttributes.pl \
  --gff "$par_gff" \
  -o "$par_output" \
  ${par_force:+--force} \
  ${par_config:+--config "${par_config}"}
