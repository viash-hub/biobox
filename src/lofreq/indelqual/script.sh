#!/bin/bash

## VIASH START
## VIASH END

# Unset all parameters that are set to "false"
[[ "$par_dindel" == "false" ]] && unset par_dindel
[[ "$par_verbose" == "false" ]] && unset par_verbose

# run lofreq indelqual
lofreq indelqual \
  -o "$par_out" \
  ${par_uniform:+-u "${par_uniform}"} \
  ${par_dindel:+--dindel} \
  ${par_ref:+-f "${par_ref}"} \
  ${par_verbose:+--verbose} \
  "$par_input"
