#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags (using loop for many parameters)
unset_if_false=(
  par_strand_windows
  par_same_strand
  par_opposite_strand
  par_unique
  par_count
  par_no_overlaps
  par_header
  par_uncompressed_bam
  par_bed_output
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Build command arguments array
cmd_args=()

# Input files (must have either input_a or input_bam)
if [[ -n "$par_input_bam" ]]; then
  cmd_args+=(-abam "$par_input_bam")
else
  cmd_args+=(-a "$par_input_a")
fi

cmd_args+=(
  -b "$par_input_b"
  ${par_window_size:+-w "$par_window_size"}
  ${par_left_window:+-l "$par_left_window"}
  ${par_right_window:+-r "$par_right_window"}
  ${par_strand_windows:+-sw}
  ${par_same_strand:+-sm}
  ${par_opposite_strand:+-Sm}
  ${par_unique:+-u}
  ${par_count:+-c}
  ${par_no_overlaps:+-v}
  ${par_header:+-header}
  ${par_uncompressed_bam:+-ubam}
  ${par_bed_output:+-bed}
)

# Execute bedtools window
bedtools window "${cmd_args[@]}" > "$par_output"
