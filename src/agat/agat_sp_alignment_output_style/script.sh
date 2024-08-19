#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# unset flags
[[ "$par_verbose" == "false" ]] && unset par_verbose

# run agat_sp_alignment_output_style.pl
agat_sp_alignment_output_style.pl \
  --gff "$par_gff" \
  --output "$par_output" \
  ${par_verbose:+-v} \
  ${par_config:+--config "${par_config}"}
