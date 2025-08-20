#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_percentage" == "false" ]] && unset par_percentage
[[ "$par_header" == "false" ]] && unset par_header

# Validate parameter combinations
# Either use -s alone, or use -p and -m together
if [[ -n "$par_shift" ]] && ([[ -n "$par_plus_shift" ]] || [[ -n "$par_minus_shift" ]]); then
  echo "ERROR: Cannot use --shift (-s) together with --plus_shift (-p) or --minus_shift (-m)" >&2
  exit 1
fi

if [[ -n "$par_plus_shift" ]] && [[ -z "$par_minus_shift" ]]; then
  echo "ERROR: --plus_shift (-p) requires --minus_shift (-m) to be specified" >&2
  exit 1
fi

if [[ -n "$par_minus_shift" ]] && [[ -z "$par_plus_shift" ]]; then
  echo "ERROR: --minus_shift (-m) requires --plus_shift (-p) to be specified" >&2
  exit 1
fi

if [[ -z "$par_shift" ]] && [[ -z "$par_plus_shift" ]] && [[ -z "$par_minus_shift" ]]; then
  echo "ERROR: Must specify either --shift (-s) or both --plus_shift (-p) and --minus_shift (-m)" >&2
  exit 1
fi

# Build command arguments array
cmd_args=(
  -i "$par_input"
  -g "$par_genome"
  ${par_shift:+-s "$par_shift"}
  ${par_plus_shift:+-p "$par_plus_shift"}
  ${par_minus_shift:+-m "$par_minus_shift"}
  ${par_percentage:+-pct}
  ${par_header:+-header}
)

# Execute bedtools shift
bedtools shift "${cmd_args[@]}" > "$par_output"
