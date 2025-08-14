#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Build the command
cmd_args=(
    # Options
    ${par_seed:+-s "$par_seed"}
    ${par_record_count:+-n "$par_record_count"}
    ${par_read_length:+--read-length "$par_read_length"}
    
    # Output files
    "$par_r1_dst"
    "$par_r2_dst"
)

# Run fq generate
fq generate "${cmd_args[@]}"
