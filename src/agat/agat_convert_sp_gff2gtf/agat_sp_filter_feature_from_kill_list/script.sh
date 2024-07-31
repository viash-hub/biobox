#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_verbose" == "false" ]] && unset par_verbose

# run agat_sp_filter_feature_from_kill_list
agat_sp_filter_feature_from_kill_list.pl \
  --gff "$par_gff" \
  --output "$par_output" \
  --kl "$par_kl" \
  ${par_type:+--type "${par_type}"} \
  ${par_attribute:+--attribute "${par_attribute}"} \
  ${par_config:+--config "${par_config}"} \
  ${par_verbose:+--verbose}
