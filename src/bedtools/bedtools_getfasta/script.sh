#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
unset_if_false=(
  par_rna
  par_strandedness
  par_split
  par_full_header
  par_name
  par_name_only
  par_tab
  par_bed_out
)

for par in "${unset_if_false[@]}"; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset "$par"
done

# Build command arguments array
cmd_args=(
    -fi "$par_input_fasta"
    -bed "$par_input_bed"
    -fo "$par_output"
    ${par_rna:+-rna}
    ${par_strandedness:+-s}
    ${par_split:+-split}
    ${par_full_header:+-fullHeader}
    ${par_name:+-name}
    ${par_name_only:+-nameOnly}
    ${par_tab:+-tab}
    ${par_bed_out:+-bedOut}
)

# Execute bedtools command
bedtools getfasta "${cmd_args[@]}"

