#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Build command arguments
args=(-i "$par_input" -g "$par_genome")

# Add flanking distance options (mutually exclusive groups)
if [ -n "$par_both" ]; then
    args+=(-b "$par_both")
elif [ -n "$par_left" ] && [ -n "$par_right" ]; then
    args+=(-l "$par_left" -r "$par_right")
elif [ -n "$par_left" ] || [ -n "$par_right" ]; then
    echo "Error: --left and --right must be used together" >&2
    exit 1
else
    echo "Error: Must specify either --both or both --left and --right" >&2
    exit 1
fi

# Add behavioral options
if [ "$par_strand" = "true" ]; then
    args+=(-s)
fi

if [ "$par_percent" = "true" ]; then
    args+=(-pct)
fi

# Add output options
if [ "$par_header" = "true" ]; then
    args+=(-header)
fi

# Execute bedtools flank
bedtools flank "${args[@]}" > "$par_output"
