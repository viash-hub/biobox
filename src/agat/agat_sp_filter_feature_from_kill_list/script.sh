#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_verbose" == "false" ]] && unset par_verbose

# convert par_type to comma separated list
par_type=$(echo $par_type | tr ';' ',').

# run agat_sp_filter_feature_from_kill_list
agat_sp_filter_feature_from_kill_list.pl \
  --gff "$par_gff" \
  --kill_list "$par_kill_list" \
  --output "$par_output" \
  ${par_type:+--type "${par_type}"} \
  ${par_attribute:+--attribute "${par_attribute}"} \
  ${par_config:+--config "${par_config}"} \
  ${par_verbose:+-v}
