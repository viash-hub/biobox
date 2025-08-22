#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_strand" == "false" ]] && unset par_strand

# Execute bedtools cluster
bedtools cluster \
  -i "$par_input" \
  ${par_distance:+-d "$par_distance"} \
  ${par_strand:+-s} \
  > "$par_output"
