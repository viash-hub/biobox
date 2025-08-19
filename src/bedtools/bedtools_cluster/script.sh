#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Unset variables that are false
unset_if_false=( par_strand )

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Execute bedtools cluster
bedtools cluster \
    -i "$par_input" \
    ${par_distance:+-d "$par_distance"} \
    ${par_strand:+-s} \
    > "$par_output"
