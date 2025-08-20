#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_strand_aware" == "false" ]] && unset par_strand_aware
[[ "$par_percentage" == "false" ]] && unset par_percentage
[[ "$par_header" == "false" ]] && unset par_header

# Validate parameter combinations
# Either use -b alone, or use -l and -r together
if [[ -n "$par_both" ]] && ([[ -n "$par_left" ]] || [[ -n "$par_right" ]]); then
  echo "ERROR: Cannot use --both (-b) together with --left (-l) or --right (-r)" >&2
  exit 1
fi

if [[ -n "$par_left" ]] && [[ -z "$par_right" ]]; then
  echo "ERROR: --left (-l) requires --right (-r) to be specified" >&2
  exit 1
fi

if [[ -n "$par_right" ]] && [[ -z "$par_left" ]]; then
  echo "ERROR: --right (-r) requires --left (-l) to be specified" >&2
  exit 1
fi

if [[ -z "$par_both" ]] && [[ -z "$par_left" ]] && [[ -z "$par_right" ]]; then
  echo "ERROR: Must specify either --both (-b) or both --left (-l) and --right (-r)" >&2
  exit 1
fi

# Build command arguments array
cmd_args=(
  -i "$par_input"
  -g "$par_genome"
  ${par_both:+-b "$par_both"}
  ${par_left:+-l "$par_left"}
  ${par_right:+-r "$par_right"}
  ${par_strand_aware:+-s}
  ${par_percentage:+-pct}
  ${par_header:+-header}
)

# Execute bedtools slop
bedtools slop "${cmd_args[@]}" > "$par_output"
