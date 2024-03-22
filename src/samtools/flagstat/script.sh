#!/bin/bash

## VIASH START
## VIASH END

set -e

output_file="$par_output"

samtools flagstat \
    "$par_bam" \
    > "$output_file"
    