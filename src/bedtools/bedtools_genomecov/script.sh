#!/bin/bash

## VIASH START
## VIASH END

[[ "$par_input_bam" == "false" ]] && unset par_input_bam

bedtools genomecov \
    ${par_input_bam:+-ibam "$par_input_bam"} \
    -i "$par_input" \
    -g "$par_genome" \
    > "$par_output"
    