#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_reciprocal" == "false" ]] && unset par_reciprocal
[[ "$par_same_strand" == "false" ]] && unset par_same_strand
[[ "$par_opposite_strand" == "false" ]] && unset par_opposite_strand
[[ "$par_include_duplicates" == "false" ]] && unset par_include_duplicates
[[ "$par_include_failed_qc" == "false" ]] && unset par_include_failed_qc
[[ "$par_proper_pairs_only" == "false" ]] && unset par_proper_pairs_only
[[ "$par_split" == "false" ]] && unset par_split

# Convert semicolon-separated bams to array
IFS=';' read -ra bams_array <<< "$par_bams"

# Build command arguments array
cmd_args=(
    -bams "${bams_array[@]}"
    -bed "$par_bed"
    ${par_min_overlap:+-f "$par_min_overlap"}
    ${par_reciprocal:+-r}
    ${par_same_strand:+-s}
    ${par_opposite_strand:+-S}
    ${par_min_quality:+-q "$par_min_quality"}
    ${par_include_duplicates:+-D}
    ${par_include_failed_qc:+-F}
    ${par_proper_pairs_only:+-p}
    ${par_split:+-split}
)

# Execute bedtools multicov
bedtools multicov "${cmd_args[@]}" > "$par_output"
