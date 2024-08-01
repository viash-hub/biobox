#!/bin/bash

## VIASH START
## VIASH END

# Unset parameters
[[ "$par_bed12" == "false" ]] && unset par_bed12
[[ "$par_uncompress_bam" == "false" ]] && unset par_uncompress_bam

# Execute bedtools bed to bam 
bedtools bedtobam \
    ${par_bed12:+-bed12} \
    ${par_uncompress_bam:+-ubam} \
    -i "$par_input" \
    -g "$par_genome" \
    > "$par_output"
