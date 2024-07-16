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

echo "> Run seqtk_subseq on FASTA file"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/a.1.fastq.gz" \
  --name_list "$meta_resources_dir/test_data/id.list" \
  --output "sub_sample.fq"

echo ">> Check if output exists"
if [ ! -f "sub_sample.fq" ]; then
    echo ">> sub_sample.fq does not exist"
    exit 1
fi

echo ">> Check number of lines in output"
n_lines=$(wc -l < sub_sample.fq)  
n_lines=$(echo "$n_lines" | awk '{print $1}')

if [ "$n_lines" -ne 2 ]; then
    echo ">> sub_sample.fq does not contain exactly two lines"
    exit 1
fi

echo ">> Check content in output"
result=$(sed -n '2p' sub_sample.fq)
expected=$(sed -n '2p' "$meta_resources_dir/test_data/a.1.fastq")
if [ "$result" == "$expected" ]; then
    echo "--> content are equal"
else
    echo "--> content are not equal"
fi

#########################################################################################
# test tab option
# echo "> Run seqtk_subseq with TAB option"
# "$meta_executable" \
#   --tab \
#   --input "$meta_resources_dir/test_data/input.fa" \
#   --name_list "$meta_resources_dir/test_data/list.lst" \
#   --output "sub_sampled.fa"

#########################################################################################
# test strand aware option
# echo "> Run seqtk_subseq with Strand Aware option"
# "$meta_executable" \
#   --strand_aware \
#   --input "$meta_resources_dir/test_data/input.fa" \
#   --name_list "$meta_resources_dir/test_data/list.lst" \
#   --output "sub_sampled.fa"

#########################################################################################
# test sequence line length option
# echo "> Run seqtk_subseq with line length option"
# "$meta_executable" \
#   --sequence_line_length 16 \
#   --input "$meta_resources_dir/test_data/input.fa" \
#   --name_list "$meta_resources_dir/test_data/list.lst" \
#   --output "sub_sampled.fa"


