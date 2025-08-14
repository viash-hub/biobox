#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_summary" == "false" ]] && unset par_summary
[[ "$par_names" == "false" ]] && unset par_names
[[ "$par_verbose_inspect" == "false" ]] && unset par_verbose_inspect
[[ "$par_large_index" == "false" ]] && unset par_large_index

# Build the command arguments
cmd_args=(
    "$par_index"
    ${par_summary:+-s}
    ${par_names:+-n}
    ${par_across:+-a "$par_across"}
    ${par_verbose_inspect:+-v}
    ${par_large_index:+--large-index}
)

# Run bowtie2-inspect
if [[ -n "$par_output" ]]; then
  bowtie2-inspect "${cmd_args[@]}" > "$par_output"
else
  bowtie2-inspect "${cmd_args[@]}"
fi
