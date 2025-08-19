#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
unset_if_false=(
  par_bedpe
  par_mate1
  par_bed12
  par_split
  par_splitD
  par_edit_distance
  par_cigar
)

for par in "${unset_if_false[@]}"; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset "$par"
done

# Build command arguments array
cmd_args=(
    -i "$par_input"
    ${par_bedpe:+-bedpe}
    ${par_mate1:+-mate1}
    ${par_bed12:+-bed12}
    ${par_split:+-split}
    ${par_splitD:+-splitD}
    ${par_edit_distance:+-ed}
    ${par_tag:+-tag "$par_tag"}
    ${par_color:+-color "$par_color"}
    ${par_cigar:+-cigar}
)

# Execute bedtools bamtobed
bedtools bamtobed "${cmd_args[@]}" > "$par_output"

