#!/bin/bash

## VIASH START
## VIASH END

set -e
[[ "$par_multiple" == "false" ]] && unset par_multiple
[[ "$par_bai" == "false" ]] && unset par_bai
[[ "$par_csi" == "false" ]] && unset par_csi
[[ "$par_multiple" == "true" ]] && par_multiple="--multiple"

samtools index \
    "$par_input" \
    ${par_csi:+-c} \
    ${par_bai:+-b} \
    ${par_min_shift:+-m "par_output_bai"} \
    ${par_multiple:+-M} \
    -o "$par_output"