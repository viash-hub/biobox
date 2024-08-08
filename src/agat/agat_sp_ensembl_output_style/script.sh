#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_verbose" == "false" ]] && unset par_verbose

# run agat_convert_sp_bed2gff.pl
agat_convert_bed2gff.pl \
  --gff "$par_bed" \
  --output "$par_output" \
  ${par_verbose:+-v}
  ${par_config:+--config "${par_config}"} \
