#!/bin/bash

## VIASH START
## VIASH END


bedtools groupby \
    -i "$par_input" \
    > "$par_output"
    