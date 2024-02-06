#!/bin/bash

## VIASH START
## VIASH END

"$meta_executable" \
    --output datasets.txt 

echo ">> Checking output"
[ ! -f "datasets.txt" ] && echo "datasets.txt does not exist" && exit 1

echo ">> Checking if output is empty"
[ ! -s "datasets.txt" ] && echo "datasets.txt is empty" && exit 1
