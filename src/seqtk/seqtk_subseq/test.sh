#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/seqtk/seqtk_subseq"
meta_resources_dir="src/seqtk"
## VIASH END

#########################################################################################
mkdir test1
cd test1

echo "> Run seqtk_subseq on FASTA/Q file"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/a.1.fastq" \
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

########################################################################################
#-- tab option --
cd ..
mkdir test2
cd test2

echo "> Run seqtk_subseq with TAB option"
"$meta_executable" \
  --tab \
  --input "$meta_resources_dir/test_data/input.fasta" \
  --name_list "$meta_resources_dir/test_data/ids.txt" \
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

cat sub_sample.fq

########################################################################################
-- strand aware option --
cd ..
mkdir test3
cd test3
echo "> Run seqtk_subseq with Strand Aware option"

"$meta_executable" \
  --strand_aware \
  --input "$meta_resources_dir/test_data/a.1.fastq" \
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

########################################################################################
-- sequence line length option --
cd ..
mkdir test4
cd test4

echo "> Run seqtk_subseq with line length option"
"$meta_executable" \
  --sequence_line_length 10 \
  --input "$meta_resources_dir/test_data/a.1.fastq" \
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
