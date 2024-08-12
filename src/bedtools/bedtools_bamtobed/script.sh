#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset parameters
[[ "$par_bedpe" == "false" ]] && unset par_bedpe
[[ "$par_mate1" == "false" ]] && unset par_mate1
[[ "$par_bed12" == "false" ]] && unset par_bed12
[[ "$par_split" == "false" ]] && unset par_split
[[ "$par_splitD" == "false" ]] && unset par_splitD
[[ "$par_edit_distance" == "false" ]] && unset par_edit_distance
[[ "$par_tag" == "false" ]] && unset par_tag
[[ "$par_color" == "false" ]] && unset par_color
[[ "$par_cigar" == "false" ]] && unset par_cigar

# Execute bedtools sort with the provided arguments
bedtools bamtobed \
    ${par_bedpe:+-bedpe} \
    ${par_mate1:+-mate1} \
    ${par_bed12:+-bed12} \
    ${par_split:+-split} \
    ${par_splitD:+-splitD} \
    ${par_edit_distance:+-ed} \
    ${par_tag:+-tag "$par_tag"} \
    ${par_cigar:+-cigar} \
    ${par_color:+-color "$par_color"} \
    -i "$par_input" \
    > "$par_output"

