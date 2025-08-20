#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_strand" == "false" ]] && unset par_strand
[[ "$par_percent" == "false" ]] && unset par_percent
[[ "$par_header" == "false" ]] && unset par_header

# Validate flanking distance options (mutually exclusive groups)
if [ -n "$par_both" ]; then
  flanking_args=(-b "$par_both")
elif [ -n "$par_left" ] && [ -n "$par_right" ]; then
  flanking_args=(-l "$par_left" -r "$par_right")
elif [ -n "$par_left" ] || [ -n "$par_right" ]; then
  echo "Error: --left and --right must be used together" >&2
  exit 1
else
  echo "Error: Must specify either --both or both --left and --right" >&2
  exit 1
fi

# Build command arguments array
cmd_args=(
  -i "$par_input"
  -g "$par_genome"
  "${flanking_args[@]}"
  ${par_strand:+-s}
  ${par_percent:+-pct}
  ${par_header:+-header}
)

# Execute bedtools flank
bedtools flank "${cmd_args[@]}" > "$par_output"
