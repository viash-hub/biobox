#!/bin/bash

## VIASH START
## VIASH END

# Unset parameters


# Execute bedtools merge with the provided arguments
bedtools merge \
    -i "$par_input" \
    > "$par_output"
