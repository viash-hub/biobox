#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Check that at least one filtering option is provided
if [ -z "$par_names" ] && [ -z "$par_sequence_pattern" ]; then
  echo "Error: Please specify either --names or --sequence_pattern for filtering." >&2
  exit 1
fi

fq filter \
    --dsts "$par_output" \
    ${par_names:+--names "${par_names}"} \
    ${par_sequence_pattern:+--sequence-pattern "${par_sequence_pattern}"} \
    "$par_input"
