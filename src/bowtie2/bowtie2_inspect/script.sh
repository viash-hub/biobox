#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_summary" == "false" ]] && unset par_summary
[[ "$par_names" == "false" ]] && unset par_names
[[ "$par_verbose_inspect" == "false" ]] && unset par_verbose_inspect
[[ "$par_debug" == "false" ]] && unset par_debug
[[ "$par_sanitized" == "false" ]] && unset par_sanitized
[[ "$par_verbose" == "false" ]] && unset par_verbose
[[ "$par_large_index" == "false" ]] && unset par_large_index

index_path="${par_index_root%/}/$par_index_prefix"

# Build the command arguments
cmd_args=(
    "$index_path"
    ${par_summary:+-s}
    ${par_names:+-n}
    ${par_across:+-a "$par_across"}
    ${par_verbose_inspect:+-v}
    ${par_debug:+--debug}
    ${par_sanitized:+--sanitized}
    ${par_verbose:+--verbose}
    ${par_large_index:+--large-index}
)

# Run bowtie2-inspect
if [[ -n "$par_output" ]]; then
  bowtie2-inspect "${cmd_args[@]}" > "$par_output"
else
  bowtie2-inspect "${cmd_args[@]}"
fi
