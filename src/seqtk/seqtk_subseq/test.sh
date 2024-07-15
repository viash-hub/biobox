#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/seqtk/seqtk_subseq"
meta_resources_dir="src/seqtk"
## VIASH END

#########################################################################################
mkdir seqtk_subseq_test
cd seqtk_subseq_test

echo "> Run seqtk_subseq on FASTA/Q file"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/input.fa" \
  --name_list "$meta_resources_dir/test_data/list.lst" \
  --output "sub_sampled.fa"

echo ">> Check if output exists"
if [ ! -f "sub_sampled.fa" ]; then
    echo ">> sub_sampled.fa does not exist"
    exit 1
fi

echo ">> Count number of subsamples in output"
num_samples=$(grep -c '^@' sub_sampled.fa)
if [ "$num_samples" -ne 2 ]; then
    echo ">> sub_sampled.fa does not contain the 2 sub-samples"
    exit 1
fi

#########################################################################################
# ... add more tests here ...

