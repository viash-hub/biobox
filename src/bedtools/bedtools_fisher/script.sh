#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Build command arguments
args=(-a "$par_input_a" -b "$par_input_b" -g "$par_genome")

# Add overlap options
if [ "$par_merge_overlaps" = "true" ]; then
    args+=(-m)
fi

if [ -n "$par_min_overlap_a" ]; then
    args+=(-f "$par_min_overlap_a")
fi

if [ -n "$par_min_overlap_b" ]; then
    args+=(-F "$par_min_overlap_b")
fi

if [ "$par_reciprocal" = "true" ]; then
    args+=(-r)
fi

if [ "$par_either" = "true" ]; then
    args+=(-e)
fi

# Add strand options
if [ "$par_same_strand" = "true" ]; then
    args+=(-s)
fi

if [ "$par_opposite_strand" = "true" ]; then
    args+=(-S)
fi

# Add format options
if [ "$par_split" = "true" ]; then
    args+=(-split)
fi

if [ "$par_bed_output" = "true" ]; then
    args+=(-bed)
fi

if [ "$par_header" = "true" ]; then
    args+=(-header)
fi

# Add advanced options
if [ "$par_no_name_check" = "true" ]; then
    args+=(-nonamecheck)
fi

if [ "$par_no_buffer" = "true" ]; then
    args+=(-nobuf)
fi

if [ -n "$par_io_buffer" ]; then
    args+=(-iobuf "$par_io_buffer")
fi

# Execute bedtools fisher
bedtools fisher "${args[@]}" > "$par_output"
