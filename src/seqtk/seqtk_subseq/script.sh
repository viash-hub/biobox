#!/bin/bash

## VIASH START
## VIASH END

[[ "$par_tab" == "false" ]] && unset par_tab
[[ "$par_strand_aware" == "false" ]] && unset par_strand_aware

seqtk subseq \
    ${par_tab:+-t} \
    ${par_strand_aware:+-s} \
    ${par_sequence_line_length:+-l "$par_sequence_line_length"} \
    "$par_input" \
    "$par_name_list" \
    > "$par_output"
