#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset parameters
[[ "$par_n_score" == "false" ]] && unset par_n_score

# Execute bedtools bed12tobed6 conversion 
bedtools bed12tobed6 \
    ${par_n_score:+-n} \
    -i "$par_input" \
    > "$par_output"
