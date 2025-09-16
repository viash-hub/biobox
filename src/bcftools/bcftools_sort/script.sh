#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Build command array
cmd_args=(
  bcftools sort
  ${par_max_mem:+--max-mem "$par_max_mem"}
  ${par_output_type:+--output-type "$par_output_type"}
  ${par_temp_dir:+--temp-dir "$par_temp_dir"}
  ${par_verbosity:+--verbosity "$par_verbosity"}
  ${par_write_index:+--write-index="$par_write_index"}
  ${par_output:+--output "$par_output"}
  "$par_input"
)

# Execute command
"${cmd_args[@]}"

