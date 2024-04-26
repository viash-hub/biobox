#!/bin/bash

set -eo pipefail

[[ "$par_nogroup" == "false" ]] && unset par_nogroup
[[ "$par_bisulfite" == "false" ]] && unset par_bisulfite
[[ "$par_reverse_compliment" == "false" ]] && unset par_reverse_compliment

IFS=";" read -ra input <<< $par_input

$(which falco) \
  ${par_nogroup:+--nogroup} \
  ${par_contaminants:+--contaminants "$par_contaminants"} \
  ${par_adapters:+--adapters "$par_adapters"} \
  ${par_limits:+--limits "$par_limits"} \
  ${par_subsample:+-subsample $par_subsample} \
  ${par_bisulfite:+-bisulfite} \
  ${par_reverse_compliment:+-reverse-compliment} \
  ${par_outdir:+--outdir "$par_outdir"} \
  ${par_format:+--format "$par_format"} \
  ${par_data_filename:+-data-filename "$par_data_filename"} \
  ${par_report_filename:+-report-filename "$par_report_filename"} \
  ${par_summary_filename:+-summary-filename "$par_summary_filename"} \
  ${input[*]}
