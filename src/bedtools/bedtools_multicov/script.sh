#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags (using loop for many parameters)
unset_if_false=(
  par_reciprocal
  par_same_strand
  par_opposite_strand
  par_include_duplicates
  par_include_failed_qc
  par_proper_pairs_only
  par_split
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Convert semicolon-separated bams to array
IFS=';' read -ra bams_array <<< "$par_bams"

# Build command arguments array
cmd_args=(
  -bams "${bams_array[@]}"
  -bed "$par_bed"
  ${par_min_overlap:+-f "$par_min_overlap"}
  ${par_reciprocal:+-r}
  ${par_same_strand:+-s}
  ${par_opposite_strand:+-S}
  ${par_min_quality:+-q "$par_min_quality"}
  ${par_include_duplicates:+-D}
  ${par_include_failed_qc:+-F}
  ${par_proper_pairs_only:+-p}
  ${par_split:+-split}
)

# Execute bedtools multicov
bedtools multicov "${cmd_args[@]}" > "$par_output"
