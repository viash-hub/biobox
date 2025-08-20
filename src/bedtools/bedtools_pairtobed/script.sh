#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_uncompressed_bam" == "false" ]] && unset par_uncompressed_bam
[[ "$par_bedpe_output" == "false" ]] && unset par_bedpe_output
[[ "$par_same_strand" == "false" ]] && unset par_same_strand
[[ "$par_opposite_strand" == "false" ]] && unset par_opposite_strand
[[ "$par_edit_distance" == "false" ]] && unset par_edit_distance

# Build command arguments array
cmd_args=()

# Handle input type - either BEDPE or BAM
if [[ -n "$par_bam_input" ]]; then
    cmd_args+=(-abam "$par_bam_input")
else
    cmd_args+=(-a "$par_bedpe")
fi

# Add BED file
cmd_args+=(-b "$par_bed")

# Add optional parameters
cmd_args+=(
    ${par_min_overlap:+-f "$par_min_overlap"}
    ${par_type:+-type "$par_type"}
    ${par_same_strand:+-s}
    ${par_opposite_strand:+-S}
    ${par_uncompressed_bam:+-ubam}
    ${par_bedpe_output:+-bedpe}
    ${par_edit_distance:+-ed}
)

# Execute bedtools pairtobed
bedtools pairtobed "${cmd_args[@]}" > "$par_output"
