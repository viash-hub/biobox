#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# unset flags
[[ "$par_bam" == "false" ]] && unset par_bam
[[ "$par_sam" == "false" ]] && unset par_sam

# run agat_convert_sp_bed2gff.pl
agat_convert_minimap2_bam2gff.pl \
  --input "$par_input" \
  --output "$par_output" \
  ${par_bam:+--bam} \
  ${par_sam:+--sam} \
  ${par_config:+--config "${par_config}"}
