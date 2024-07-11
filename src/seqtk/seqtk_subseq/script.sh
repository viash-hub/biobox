#!/bin/bash

## VIASH START
## VIASH END

seqtk sample \
    ${par_tab:+-t} \
    ${par_strand_aware:+-s} \
    ${par_line_length:+-l "$par_line_length"} \
    "$par_input" \
    > "$par_output"