#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Create temporary directory for bcftools if meta_temp_dir is available
if [[ -n "$meta_temp_dir" ]]; then
  bcftools_temp_dir=$(mktemp -d "$meta_temp_dir/bcftools_sort.XXXXXX")
  # Set up cleanup trap
  trap 'rm -rf "$bcftools_temp_dir"' EXIT
fi

# Build command array
cmd_args=(
  bcftools sort
  ${meta_memory_mb:+--max-mem "${meta_memory_mb}M"}
  ${par_output_type:+--output-type "$par_output_type"}
  ${bcftools_temp_dir:+--temp-dir "$bcftools_temp_dir"}
  ${par_verbosity:+--verbosity "$par_verbosity"}
  ${par_write_index:+--write-index="$par_write_index"}
  ${par_output:+--output "$par_output"}
  "$par_input"
)

# Execute command
"${cmd_args[@]}"

