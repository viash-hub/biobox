#!/bin/bash

## VIASH START
## VIASH END

seqtk subseq \
    ${par_tab:+-t} \
    ${par_strand_aware:+-s} \
    ${par_sequence_line_length:+-l "$par_sequence_line_length"} \
    "$par_input" \
    "$par_name_list" \
    > "$par_output"