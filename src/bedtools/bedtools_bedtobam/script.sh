#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters
[[ "$par_bed12" == "false" ]] && unset par_bed12
[[ "$par_uncompress_bam" == "false" ]] && unset par_uncompress_bam

# Build command arguments array
cmd_args=(
  -i "$par_input"
  -g "$par_genome"
  ${par_map_quality:+-mapq "$par_map_quality"}
  ${par_bed12:+-bed12}
  ${par_uncompress_bam:+-ubam}
)

# Execute bedtools bedtobam
bedtools bedtobam "${cmd_args[@]}" > "$par_output"
