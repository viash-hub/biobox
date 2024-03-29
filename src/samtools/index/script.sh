#!/bin/bash

## VIASH START
## VIASH END

set -e

[[ "$par_output_csi" == "false" ]] && unset par_output_csi

samtools index \
    "$par_input" \
    ${par_output_csi:+-c} \
    -o "$par_output"