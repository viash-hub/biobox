#!/bin/bash

## VIASH START
## VIASH END

# unset flags
[[ "$par_verbose" == "false" ]] && unset par_verbose

# Convert a list of file names to multiple -gff arguments
input_files=""
IFS=";" read -ra file_names <<< "$par_gff"
for file in "${file_names[@]}"; do
    input_files+="--gff $file "
done
unset IFS


# run agat_sp_complement_annotations.pl
agat_sp_complement_annotations.pl \
  --ref "$par_ref" \
  $input_files \
  -o "$par_output" \
  ${par_size_min:+--size_min "${par_size_min}"} \
  ${par_config:+--config "${par_config}"} \
  ${par_verbose:+--verbose}
