#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# exclusive OR for required arguments $par_probability and $par_record_count
if { [ -n "$par_probability" ] && [ -n "$par_record_count" ]; } || \
   { [ -z "$par_probability" ] && [ -z "$par_record_count" ]; }; then
  echo "Error: Please specify either --probability or --record_count, but not both." >&2
  exit 1
fi

fq subsample \
    ${par_output_1:+--r1-dst "${par_output_1}"} \
    ${par_output_2:+--r2-dst "${par_output_2}"} \
    ${par_probability:+--probability "${par_probability}"} \
    ${par_record_count:+--record-count "${par_record_count}"} \
    ${par_seed:+--seed "${par_seed}"} \
    ${par_input_1} \
    ${par_input_2:+"$par_input_2"}
