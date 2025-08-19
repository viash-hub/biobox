#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Unset variables that are false
unset_if_false=( par_limit_chromosomes )

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

# Execute bedtools complement
bedtools complement \
    -i "$par_input" \
    -g "$par_genome" \
    ${par_limit_chromosomes:+-L} \
    > "$par_output"
