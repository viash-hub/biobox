#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Build the command
cmd_args=(
    # Options
    ${par_max_occ:+-n "$par_max_occ"}
    ${par_output:+-f "$par_output"}
    ${par_read_group:+-r "$par_read_group"}
    
    # Required arguments: index, SAI file, FASTQ file
    "$par_index"
    "$par_sai"
    "$par_reads"
)

# Run bwa samse
bwa samse "${cmd_args[@]}"
