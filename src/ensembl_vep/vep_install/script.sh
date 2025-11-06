#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters (biobox standard)
[[ "$par_no_update" == "false" ]] && unset par_no_update
[[ "$par_no_plugins" == "false" ]] && unset par_no_plugins
[[ "$par_no_htslib" == "false" ]] && unset par_no_htslib
[[ "$par_no_bioperl" == "false" ]] && unset par_no_bioperl
[[ "$par_prefer_bin" == "false" ]] && unset par_prefer_bin
[[ "$par_convert" == "false" ]] && unset par_convert
[[ "$par_quiet" == "false" ]] && unset par_quiet

# Handle multiple values (semicolon-separated from Viash)
if [[ -n "$par_species" ]]; then
  # Convert semicolon-separated to comma-separated for VEP
  par_species=$(echo "$par_species" | tr ';' ',')
fi

if [[ -n "$par_plugins" ]]; then
  # Convert semicolon-separated to comma-separated for VEP
  par_plugins=$(echo "$par_plugins" | tr ';' ',')
fi

# Build command array (preferred pattern)
cmd_args=(
  vep_install
  ${par_destdir:+--destdir "$par_destdir"}
  ${par_cachedir:+--cachedir "$par_cachedir"}
  ${par_cache_version:+--cache_version "$par_cache_version"}
  ${par_auto:+--auto "$par_auto"}
  ${par_species:+--species "$par_species"}
  ${par_assembly:+--assembly "$par_assembly"}
  ${par_plugins:+--plugins "$par_plugins"}
  ${par_pluginsdir:+--pluginsdir "$par_pluginsdir"}
  ${par_no_update:+--no_update}
  ${par_no_plugins:+--no_plugins}
  ${par_no_htslib:+--no_htslib}
  ${par_no_bioperl:+--no_bioperl}
  ${par_prefer_bin:+--prefer_bin}
  ${par_convert:+--convert}
  ${par_quiet:+--quiet}
)

# Execute command
"${cmd_args[@]}"