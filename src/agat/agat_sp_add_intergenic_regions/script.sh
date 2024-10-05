#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# unset flags
[[ "$par_verbose" == "false" ]] && unset par_verbose

# run agat_sp_add_intergenic_regions.pl
agat_sp_add_intergenic_regions.pl \
  --gff "$par_gff" \
  --output "$par_output" \
  ${par_verbose:+--verbose} \
  ${par_config:+--config "${par_config}"}
