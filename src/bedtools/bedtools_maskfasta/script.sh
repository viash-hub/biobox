#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_soft_mask" == "false" ]] && unset par_soft_mask
[[ "$par_full_header" == "false" ]] && unset par_full_header

# Build command arguments array
cmd_args=(
    -fi "$par_input_fasta"
    -bed "$par_input_bed"
    -fo "$par_output"
    ${par_mask_character:+-mc "$par_mask_character"}
    ${par_soft_mask:+-soft}
    ${par_full_header:+-fullHeader}
)

# Execute bedtools maskfasta
bedtools maskfasta "${cmd_args[@]}"
