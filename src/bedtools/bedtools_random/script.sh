#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Build command arguments array
cmd_args=(
  -g "$par_genome"
  ${par_length:+-l "$par_length"}
  ${par_number:+-n "$par_number"}
  ${par_seed:+-seed "$par_seed"}
)

# Execute bedtools random
bedtools random "${cmd_args[@]}" > "$par_output"
