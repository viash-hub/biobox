#!/bin/bash

## VIASH START
## VIASH END

# Convert a list of file names to multiple -gff arguments
input_files=""
IFS=";" read -ra file_names <<< "$par_gff"
for file in "${file_names[@]}"; do
    input_files+="--gff $file "
done

# run agat_sp_merge_annotations
agat_sp_merge_annotations.pl \
  $input_files \
  -o "$par_output" \
  ${par_config:+--config "${par_config}"}
