#!/bin/bash

## VIASH START
## VIASH END

set -e

samtools flagstat \
    "$par_bam" \
    > "$par_output"
    