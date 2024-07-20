#!/bin/bash

test_dir="${meta_resources_dir}/test_data"

echo ">>> Testing $meta_functionality_name"

############################################################################################################

echo ">>> Test 1: Testing for paired-end reads"
"$meta_executable" \
    --input "$test_dir/scrb_seq_fastq.1_30"\
    --read2_in "$test_dir/scrb_seq_fastq.2_30" \
    --bc_pattern "CCCCCCNNNNNNNNNN"\
    --bc_pattern2 "CCCCCCNNNNNNNNNN" \
    --extract_method string \
    --umi_separator '_' \
    --grouping_method directional \
    --umi_discard_read 0 \
    --output scrb_seq_fastq.1_30.extract \
    --read2_out scrb_seq_fastq.2_30.extract \
    --random_seed 1

echo ">> Checking if the correct files are present"
[[ ! -f "scrb_seq_fastq.1_30.extract" ]] || [[ ! -f "scrb_seq_fastq.2_30.extract" ]] && echo "Reads file missing" && exit 1
[ ! -s "scrb_seq_fastq.1_30.extract" ] && echo "Read 1 file is empty" && exit 1
[ ! -s "scrb_seq_fastq.2_30.extract" ] && echo "Read 2 file is empty" && exit 1


echo ">> Checking if the files are correct"
diff -q "${meta_resources_dir}/scrb_seq_fastq.1_30.extract" "$test_dir/scrb_seq_fastq.1_30.extract" || \
    (echo "Read 1 file is not correct" && exit 1)
diff -q "${meta_resources_dir}/scrb_seq_fastq.2_30.extract" "$test_dir/scrb_seq_fastq.2_30.extract" || \
    (echo "Read 2 file is not correct" && exit 1)

rm scrb_seq_fastq.1_30.extract scrb_seq_fastq.2_30.extract

############################################################################################################

echo ">>> Test 2: Testing for paired-end reads with umi_discard_reads option"
"$meta_executable" \
    --input "$test_dir/scrb_seq_fastq.1_30" \
    --read2_in "$test_dir/scrb_seq_fastq.2_30" \
    --bc_pattern CCCCCCNNNNNNNNNN \
    --bc_pattern2 CCCCCCNNNNNNNNNN \
    --extract_method string \
    --umi_separator '_' \
    --grouping_method directional \
    --umi_discard_read 2 \
    --output scrb_seq_fastq.1_30.extract \
    --random_seed 1

echo ">> Checking if the correct files are present"
[ ! -f "scrb_seq_fastq.1_30.extract" ] && echo "Read 1 file is missing" && exit 1
[ ! -s "scrb_seq_fastq.1_30.extract" ] && echo "Read 1 file is empty" && exit 1
[ -f "scrb_seq_fastq.2_30.extract" ] && echo "Read 2 is not discarded" && exit 1

echo ">> Checking if the files are correct"
diff -q "${meta_resources_dir}/scrb_seq_fastq.1_30.extract" "$test_dir/scrb_seq_fastq.1_30.extract" || \
    (echo "Read 1 file is not correct" && exit 1)

rm scrb_seq_fastq.1_30.extract

############################################################################################################

echo ">>> Test 3: Testing for single-end reads"
"$meta_executable" \
    --input "$test_dir/slim_30.fastq" \
    --bc_pattern "^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})" \
    --extract_method regex \
    --umi_separator '_' \
    --grouping_method directional \
    --output slim_30.extract \
    --random_seed 1 

echo ">> Checking if the correct files are present"
[ ! -f "slim_30.extract" ] && echo "Trimmed reads file missing" && exit 1
[ ! -s "slim_30.extract" ] && echo "Trimmed reads file is empty" && exit 1

echo ">> Checking if the files are correct"
diff -q "${meta_resources_dir}/slim_30.extract" "$test_dir/slim_30.extract" || \
    (echo "Trimmed reads file is not correct" && exit 1)

rm slim_30.extract

echo ">>> Test finished successfully"
exit 0