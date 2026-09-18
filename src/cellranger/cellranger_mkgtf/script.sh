#!/bin/bash

set -eo pipefail

## VIASH START
par_input_gtf="test_data/reference_small.gtf"
par_output_gtf="reference_small_filtered.gtf"
par_attribute="gene_biotype:protein_coding"
## VIASH END

cmd_args=()

# `--attribute` can be specified multiple times, so split the `;`-separated
# values into one `--attribute` argument each
if [[ -n "$par_attribute" ]]; then
  IFS=';' read -ra attributes <<< "$par_attribute"
  for attribute in "${attributes[@]}"; do
    cmd_args+=("--attribute=$attribute")
  done
fi

echo "> Running cellranger mkgtf"
cellranger mkgtf \
  "$par_input_gtf" \
  "$par_output_gtf" \
  "${cmd_args[@]}"
