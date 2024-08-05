#!/bin/bash

## VIASH START
par_input='./src/nanoplot/test_data/real.fastq'
## VIASH END

NanoPlot \
    "$par_input" \
    > "$par_output"