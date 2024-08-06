#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_verbose" == "false" ]] && unset par_verbose

# run agat_sp_add_intergenic_regions.pl
agat_sp_add_intergenic_regions.pl \
  --bed "$par_bed" \
  -o "$par_output" \
  ${par_verbose:+--verbose} \
  ${par_config:+--config "${par_config}"}
