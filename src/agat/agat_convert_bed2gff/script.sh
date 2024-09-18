#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_inflate_off" == "true" ]] && unset par_inflate_off
[[ "$par_verbose" == "false" ]] && unset par_verbose

# run agat_convert_sp_bed2gff.pl
agat_convert_bed2gff.pl \
  --bed "$par_bed" \
  -o "$par_output" \
  ${par_source:+--source "${par_source}"} \
  ${par_primary_tag:+--primary_tag "${par_primary_tag}"} \
  ${par_inflate_off:+--inflate_off} \
  ${par_inflate_type:+--inflate_type "${par_inflate_type}"} \
  ${par_verbose:+--verbose}
  ${par_config:+--config "${par_config}"} \
