#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_strand_slop" == "false" ]] && unset par_strand_slop
[[ "$par_ignore_strand" == "false" ]] && unset par_ignore_strand
[[ "$par_require_different_names" == "false" ]] && unset par_require_different_names

# Build command arguments array
cmd_args=(
    -a "$par_bedpe_a"
    -b "$par_bedpe_b"
    ${par_min_overlap:+-f "$par_min_overlap"}
    ${par_type:+-type "$par_type"}
    ${par_slop:+-slop "$par_slop"}
    ${par_strand_slop:+-ss}
    ${par_ignore_strand:+-is}
    ${par_require_different_names:+-rdn}
)

# Execute bedtools pairtopair
bedtools pairtopair "${cmd_args[@]}" > "$par_output"
