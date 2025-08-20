#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags (using loop for many parameters)
unset_if_false=(
  par_sizeA
  par_sizeD
  par_chrThenSizeA
  par_chrThenSizeD
  par_chrThenScoreA
  par_chrThenScoreD
  par_header
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Execute bedtools sort
bedtools sort \
  -i "$par_input" \
  ${par_sizeA:+-sizeA} \
  ${par_sizeD:+-sizeD} \
  ${par_chrThenSizeA:+-chrThenSizeA} \
  ${par_chrThenSizeD:+-chrThenSizeD} \
  ${par_chrThenScoreA:+-chrThenScoreA} \
  ${par_chrThenScoreD:+-chrThenScoreD} \
  ${par_genome:+-g "$par_genome"} \
  ${par_faidx:+-faidx "$par_faidx"} \
  ${par_header:+-header} \
  > "$par_output"
