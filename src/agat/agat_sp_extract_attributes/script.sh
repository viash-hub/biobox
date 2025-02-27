#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# unset flags
[[ "$par_merge" == "false" ]] && unset par_merge
[[ "$par_dot" == "false" ]] && unset par_dot

# run agat_sp_extract_attributes.pl
agat_sp_extract_attributes.pl \
  --gff "$par_gff" \
  --attribute "$par_attribute" \
  --output "$par_output" \
  ${par_primary_tag:+-p "${par_primary_tag}"} \
  ${par_merge:+--merge} \
  ${par_dot:+-d} \
  ${par_config:+--config "${par_config}"}