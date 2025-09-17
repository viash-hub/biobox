#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters
[[ "$par_n_score" == "false" ]] && unset par_n_score

# Build command arguments array
cmd_args=(
  -i "$par_input"
  ${par_n_score:+-n}
)

# Execute bedtools bed12tobed6
bedtools bed12tobed6 "${cmd_args[@]}" > "$par_output"
