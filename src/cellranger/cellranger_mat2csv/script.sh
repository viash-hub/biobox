#!/bin/bash

set -eo pipefail

## VIASH START
par_input="filtered_feature_bc_matrix.h5"
par_output="matrix.csv"
## VIASH END

# cellranger mat2csv does not create the directory holding the output file itself.
mkdir -p "$(dirname "$par_output")"

cmd_args=(
  "$par_input"
  "$par_output"
  ${par_genome:+--genome="$par_genome"}
)

echo "> Running cellranger mat2csv"
cellranger mat2csv "${cmd_args[@]}"

exit 0
