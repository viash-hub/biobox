#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# unset flags
[[ "$par_d" == "false" ]] && unset par_d

if [[ -n "$par_gs_size" && -n "$par_gs_fasta" ]]; then
  echo "[error] Please provide only one of the following options to set genome size: --gs_size or --gs_fasta"
  exit 1
fi

# run agat_sp_statistics
agat_sp_statistics.pl \
  -i "$par_gff" \
  -o "$par_output" \
  ${par_plot:+-d} \
  ${par_gs_size:+--gs "${par_gs_size}"} \
  ${par_gs_fasta:+--gs "${par_gs_fasta}"} \
  ${par_verbose:+--verbose "${par_verbose}"} \
  ${par_config:+--config "${par_config}"}


