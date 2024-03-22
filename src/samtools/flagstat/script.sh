#!/bin/bash

## VIASH START
## VIASH END

set -e

samtools flagstat \
    "$par_input_bam" \
    > "$par_output_stats"
    