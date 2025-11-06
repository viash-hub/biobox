#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters (biobox standard)
[[ "$par_only_matched" == "false" ]] && unset par_only_matched
[[ "$par_ontology" == "false" ]] && unset par_ontology
[[ "$par_count" == "false" ]] && unset par_count
[[ "$par_force_valid_header" == "false" ]] && unset par_force_valid_header

# Handle multiple filter values (semicolon-separated from Viash)
if [[ -n "$par_filter" ]]; then
  IFS=';' read -ra filter_array <<< "$par_filter"
fi

# Build command array (preferred pattern)
cmd_args=(
  filter_vep
  --input_file "$par_input_file"
  --output_file "$par_output_file"
  ${par_format:+--format "$par_format"}
  ${par_only_matched:+--only_matched}
  ${par_ontology:+--ontology}
  ${par_count:+--count}
  ${par_force_valid_header:+--force_valid_header}
  ${par_test:+--test "$par_test"}
)

# Handle multiple filter conditions
if [[ -n "$par_filter" ]]; then
  for filter_condition in "${filter_array[@]}"; do
    cmd_args+=(--filter "$filter_condition")
  done
fi

# Execute command
"${cmd_args[@]}"