#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_limit_chromosomes" == "false" ]] && unset par_limit_chromosomes

# Execute bedtools complement
bedtools complement \
  -i "$par_input" \
  -g "$par_genome" \
  ${par_limit_chromosomes:+-L} \
  > "$par_output"
