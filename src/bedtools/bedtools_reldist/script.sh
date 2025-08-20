#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_detail" == "false" ]] && unset par_detail

# Build command arguments array
cmd_args=(
    -a "$par_bed_a"
    -b "$par_bed_b"
    ${par_detail:+-detail}
)

# Execute bedtools reldist
bedtools reldist "${cmd_args[@]}" > "$par_output"
