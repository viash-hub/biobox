#!/bin/bash

## VIASH START
## VIASH END

seqtk sample \
    ${par_two_pass_mode:+-2} \
    ${par_seed:+-s "$par_seed"} \
    "$par_input" \
    "$par_fraction_number" \
    > "$par_output"