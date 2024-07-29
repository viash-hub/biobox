#!/bin/bash

## VIASH START
## VIASH END

# Unset parameters
[[ "$par_sizeA" == "false" ]] && unset par_sizeA
[[ "$par_sizeD" == "false" ]] && unset par_sizeD
[[ "$par_chrThenSizeA" == "false" ]] && unset par_chrThenSizeA
[[ "$par_chrThenSizeD" == "false" ]] && unset par_chrThenSizeD
[[ "$par_chrThenScoreA" == "false" ]] && unset par_chrThenScoreA
[[ "$par_chrThenScoreD" == "false" ]] && unset par_chrThenScoreD
[[ "$par_header" == "false" ]] && unset par_header

# Execute bedtools sort with the provided arguments
bedtools sort \
    ${par_sizeA:+-sizeA} \
    ${par_sizeD:+-sizeD} \
    ${par_chrThenSizeA:+-chrThenSizeA} \
    ${par_chrThenSizeD:+-chrThenSizeD} \
    ${par_chrThenScoreA:+-chrThenScoreA} \
    ${par_chrThenScoreD:+-chrThenScoreD} \
    ${par_genome:+-g "$par_genome"} \
    ${par_faidx:+-faidx "$par_faidx"} \
    ${par_header:+-header} \
    -i "$par_input" \
    > "$par_output"
