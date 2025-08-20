#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_ubam" == "false" ]] && unset par_ubam

# Build command arguments array
cmd_args=(
  -i "$par_input"
  -g "$par_genome"
  ${par_mapq:+-mapq "$par_mapq"}
  ${par_ubam:+-ubam}
)

# Execute bedtools bedpetobam
bedtools bedpetobam "${cmd_args[@]}" > "$par_output"
