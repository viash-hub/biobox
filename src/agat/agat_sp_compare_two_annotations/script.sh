#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# unset flags
[[ "$par_debug" == "false" ]] && unset par_debug
[[ "$par_verbose" == "false" ]] && unset par_verbose

# Debug statement to check if par_output_dir is set
echo "par_output_dir is set to: ${par_output_dir}"

# Ensure par_output_dir is set
if [ -z "$par_output_dir" ]; then
  echo "Error: par_output_dir is not set."
  exit 1
fi

# run agat_sp_compare_two_annotations.pl
agat_sp_compare_two_annotations.pl \
  -gff1 "$par_gff1" \
  -gff2 "$par_gff2" \
  ${par_output_dir:+-o "${par_output_dir}"} \
  ${par_debug:+--debug} \
  ${par_verbose:+--verbose} \
  ${par_config:+--config "${par_config}"}
