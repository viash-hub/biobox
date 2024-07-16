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
#
# test fq file

# TODO: Figure out how the test fq file should look like and how the reg.bed file should look like

#########################################################################################
# test tab option
echo "> Run seqtk_subseq with TAB option"
"$meta_executable" \
  --tab \
  --input "$meta_resources_dir/test_data/input.fa" \
  --name_list "$meta_resources_dir/test_data/list.lst" \
  --output "sub_sampled.fa"

#########################################################################################
# test strand aware option
echo "> Run seqtk_subseq with Strand Aware option"
"$meta_executable" \
  --strand_aware \
  --input "$meta_resources_dir/test_data/input.fa" \
  --name_list "$meta_resources_dir/test_data/list.lst" \
  --output "sub_sampled.fa"

#########################################################################################
# test sequence line length option
echo "> Run seqtk_subseq with line length option"
"$meta_executable" \
  --sequence_line_length 16 \
  --input "$meta_resources_dir/test_data/input.fa" \
  --name_list "$meta_resources_dir/test_data/list.lst" \
  --output "sub_sampled.fa"


# echo ">> Compare reads"
# # Extract headers
# headers1=$(grep '^@' sampled_1.fastq | sed -e's/ 1$//' | sort)
# headers2=$(grep '^@' sampled_2.fastq | sed -e 's/ 2$//' | sort)

# # Compare headers
# diff <(echo "$headers1") <(echo "$headers2") || { echo "Mismatch detected"; exit 1; }