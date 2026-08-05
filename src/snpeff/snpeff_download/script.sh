#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

[[ "$par_debug" == "false" ]] && unset par_debug
[[ "$par_no_log" == "false" ]] && unset par_no_log
[[ "$par_quiet" == "false" ]] && unset par_quiet
[[ "$par_verbose" == "false" ]] && unset par_verbose

mkdir -p "$par_output"

cmd_args=(
  ${par_config:+-config "$par_config"}
  ${par_config_option:+-configOption "$par_config_option"}
  ${par_debug:+-debug}
  -dataDir "$par_output"
  ${par_no_log:+-noLog}
  ${par_quiet:+-quiet}
  ${par_verbose:+-verbose}
)

snpEff download "${cmd_args[@]}" "$par_genome_version"
