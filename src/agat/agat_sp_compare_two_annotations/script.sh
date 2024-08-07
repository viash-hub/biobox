#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_debug" == "false" ]] && unset par_debug
[[ "$par_verbose" == "false" ]] && unset par_verbose

# run agat_sp_compare_two_annotations.pl
agat_sp_compare_two_annotations.pl \
  -gff1 "$par_gff1" \
  -gff2 "$par_gff2" \
  -output "$par_output" \
  ${par_debug:+--debug} \
  ${par_verbose:+--verbose} \
  ${par_config:+--config "${par_config}"}
