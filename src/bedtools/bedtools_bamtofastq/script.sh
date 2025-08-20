#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters
[[ "$par_tags" == "false" ]] && unset par_tags

# Build command arguments array
cmd_args=(
  -i "$par_input"
  -fq "$par_fastq"
  ${par_fastq2:+-fq2 "$par_fastq2"}
  ${par_tags:+-tags}
)

# Execute bedtools bamtofastq
bedtools bamtofastq "${cmd_args[@]}"