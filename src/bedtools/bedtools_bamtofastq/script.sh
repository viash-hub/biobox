#!/bin/bash

## VIASH START
## VIASH END

# Unset parameters
[[ "$par_tags" == "false" ]] && unset par_tags

# Execute bedtools bamtofastq with the provided arguments
bedtools bamtofastq \
    ${par_tags:+-tags} \
    ${par_output_fq2:+-fq2 "$par_output_fq2"} \
    -i "$par_input" \
    -fq "$par_output_fq"
    

