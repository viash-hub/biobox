#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

if [ $par_strandedness == 'forward' ]; then
    strandedness='--strandedness forward'
elif [ $par_strandedness == 'reverse' ]; then
    strandedness="--strandedness reverse"
else
    strandedness=''
fi

IFS=";" read -ra input <<< $par_input

INDEX=$(find -L $meta_resources_dir/$par_index -name "*.grp" | sed 's/\.grp$//')

rsem-calculate-expression \
    $strandedness \
    ${par_paired:+--paired-end} \
    $par_extra_args \
    ${input[*]} \
    $INDEX \
    $par_id \
    --quiet
   
