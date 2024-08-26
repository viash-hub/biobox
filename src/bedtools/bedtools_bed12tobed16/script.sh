#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Execute bedtools bed12tobed6 conversion 
bedtools bed12tobed6 \
    -i "$par_input" \
    > "$par_output"
