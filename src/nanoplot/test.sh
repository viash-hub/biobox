#!/bin/bash

set -e

## VIASH START
meta_executable="$PWD/target/executable/NanoPlot/NanoPlot"
meta_resources_dir="$PWD/src/nanoplot"
## VIASH END

echo "> Run NanoPlot"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/real.fastq" \
  -o "output"

echo ">> Check if output directory is empty"
if [ -z "$(ls -A ./output)" ]; then
    echo "No output"
    exit 1
fi

# Check whether the contents of the resulting file is correct