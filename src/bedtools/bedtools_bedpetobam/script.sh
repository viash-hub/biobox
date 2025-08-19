#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Unset variables that are false
unset_if_false=( par_ubam )

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Execute bedtools bedpetobam
bedtools bedpetobam \
    -i "$par_input" \
    -g "$par_genome" \
    ${par_mapq:+-mapq "$par_mapq"} \
    ${par_ubam:+-ubam} \
    > "$par_output"
