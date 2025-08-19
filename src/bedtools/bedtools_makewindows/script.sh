#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_reverse" == "false" ]] && unset par_reverse

# Validate input options (mutually exclusive)
if [[ -n "$par_genome" && -n "$par_input" ]]; then
    echo "Error: Cannot use both --genome and --input options. Choose one." >&2
    exit 1
elif [[ -z "$par_genome" && -z "$par_input" ]]; then
    echo "Error: Must specify either --genome or --input option." >&2
    exit 1
fi

# Validate window options (mutually exclusive)
if [[ -n "$par_window_size" && -n "$par_num_windows" ]]; then
    echo "Error: Cannot use both --window_size and --num_windows. Choose one." >&2
    exit 1
elif [[ -z "$par_window_size" && -z "$par_num_windows" ]]; then
    echo "Error: Must specify either --window_size or --num_windows." >&2
    exit 1
fi

# Validate step size usage
if [[ -n "$par_step_size" && -z "$par_window_size" ]]; then
    echo "Error: --step_size can only be used with --window_size." >&2
    exit 1
fi

# Build command arguments array
cmd_args=(
    ${par_genome:+-g "$par_genome"}
    ${par_input:+-b "$par_input"}
    ${par_window_size:+-w "$par_window_size"}
    ${par_step_size:+-s "$par_step_size"}
    ${par_num_windows:+-n "$par_num_windows"}
    ${par_id_type:+-i "$par_id_type"}
    ${par_reverse:+-reverse}
)

# Execute bedtools makewindows
bedtools makewindows "${cmd_args[@]}" > "$par_output"
