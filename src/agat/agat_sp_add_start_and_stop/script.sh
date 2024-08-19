#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# unset flags
[[ "$par_ni" == "false" ]] && unset par_ni
[[ "$par_verbose" == "false" ]] && unset par_verbose
[[ "$par_extend" == "false" ]] && unset par_extend

# run agat_sp_add_start_and_stop.pl
agat_sp_add_start_and_stop.pl \
  --gff "$par_gff" \
  --fasta "$par_fasta" \
  --output "$par_output" \
  ${par_ct:+--ct "${par_ct}"} \
  ${par_extend:+--extend} \
  ${par_ni:+--ni} \
  ${par_verbose:+--v} \
  ${par_config:+--config "${par_config}"}
