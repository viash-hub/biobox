#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_ni" == "true" ]] && unset par_ni
[[ "$par_verbose" == "false" ]] && unset par_verbose

# run agat_convert_sp_bed2gff.pl
agat_sp_add_start_and_stop.pl \
  --gff "$par_gff" \
  --fasta "$par_fasta" \
  --output "$par_output" \
  ${par_ct:+--ct "${par_ct}"} \
  ${par_extend:+--extend} \
  ${par_ni:+--ni} \
  ${par_verbose:+--verbose} \
  ${par_config:+--config "${par_config}"}
