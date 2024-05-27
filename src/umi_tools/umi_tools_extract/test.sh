#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/out"

mkdir -p "$out_dir"

#############################################################################################

echo ">>> Test 1: Testing for paired-end reads"

"$meta_executable" \
    --input "$test_dir/scrb_seq_fastq.1.subsample.fastq" \
    --read2_in "$test_dir/scrb_seq_fastq.2.subsample.fastq" \
    --bc_pattern CCCCCCNNNNNNNNNN \
    --bc_pattern2 CCCCCCNNNNNNNNNN \
    --extract_method string \
    --output "$out_dir/scrb_seq_fastq.1.umi_extract.fastq"

echo ">> Checking if the correct files are present"
[[ ! -f "$out_dir/scrb_seq_fastq.1.umi_extract.fastq" ]] && echo "Reads file missing" && exit 1
[ ! -s "$out_dir/scrb_seq_fastq.1.umi_extract.fastq" ] && echo "Read 1 file is empty" && exit 1

rm "$out_dir/scrb_seq_fastq.1.umi_extract.fastq"

#############################################################################################

echo ">>> Test 2: Testing for paired-end reads with umi_discard_reads option"
"$meta_executable" \
    --input "$test_dir/scrb_seq_fastq.1.subsample.fastq" \
    --read2_in "$test_dir/scrb_seq_fastq.2.subsample.fastq" \
    --bc_pattern CCCCCCNNNNNNNNNN \
    --bc_pattern2 CCCCCCNNNNNNNNNN \
    --extract_method string \
    --output "$out_dir/scrb_seq_fastq.1.umi_extract.fastq"

echo ">> Checking if the correct files are present"
[ ! -f "$out_dir/scrb_seq_fastq.1.umi_extract.fastq" ] && echo "Read 1 file is missing" && exit 1
[ ! -s "$out_dir/scrb_seq_fastq.1.umi_extract.fastq" ] && echo "Read 1 file is empty" && exit 1

#############################################################################################

echo ">>> Test 3: Testing for single-end reads"
"$meta_executable" \
    --input "$test_dir/slim.subsample.fastq" \
    --bc_pattern  "^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})" \
    --extract_method regex \
    --output "$out_dir/slim.umi_extract.fastq"

echo ">> Checking if the correct files are present"
[ ! -f "$out_dir/slim.umi_extract.fastq" ] && echo "Trimmed reads file missing" && exit 1
[ ! -s "$out_dir/slim.umi_extract.fastq" ] && echo "Trimmed reads file is empty" && exit 1


echo ">>> Test finished successfully"
exit 0