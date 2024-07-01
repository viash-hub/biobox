#!/bin/bash

test_dir="${meta_resources_dir}/test_data"

echo ">>> Testing $meta_functionality_name"

############################################################################################################

echo ">>> Test 1: Testing for paired-end reads"
"$meta_executable" \
    --paired \
    --input "$test_dir/scrb_seq_fastq.1.gz"\
    --read2_in "$test_dir/scrb_seq_fastq.2.gz" \
    --bc_pattern "CCCCCCNNNNNNNNNN"\
    --bc_pattern2 "CCCCCCNNNNNNNNNN" \
    --umitools_extract_method string \
    --umitools_umi_separator '_' \
    --umitools_grouping_method directional \
    --umi_discard_read 0 \
    --read1_out scrb_seq_fastq.1.umi_extract.fastq.gz \
    --read2_out scrb_seq_fastq.2.umi_extract.fastq.gz \
    --random_seed 1

echo ">> Checking if the correct files are present"
[[ ! -f "scrb_seq_fastq.1.umi_extract.fastq.gz" ]] || [[ ! -f "scrb_seq_fastq.2.umi_extract.fastq.gz" ]] && echo "Reads file missing" && exit 1
[ ! -s "scrb_seq_fastq.1.umi_extract.fastq.gz" ] && echo "Read 1 file is empty" && exit 1
[ ! -s "scrb_seq_fastq.2.umi_extract.fastq.gz" ] && echo "Read 2 file is empty" && exit 1


echo ">> Checking if the files are correct"
diff -q <(gunzip -c "${meta_resources_dir}/scrb_seq_fastq.1.umi_extract.fastq.gz") <(gunzip -c "$test_dir/scrb_seq_fastq.1.umi_extract.fastq.gz") || \
    (echo "Read 1 file is not correct" && exit 1)
diff -q <(gunzip -c "${meta_resources_dir}/scrb_seq_fastq.2.umi_extract.fastq.gz") <(gunzip -c "$test_dir/scrb_seq_fastq.2.umi_extract.fastq.gz") || \
    (echo "Read 2 file is not correct" && exit 1)

rm scrb_seq_fastq.1.umi_extract.fastq.gz scrb_seq_fastq.2.umi_extract.fastq.gz

############################################################################################################

echo ">>> Test 2: Testing for paired-end reads with umi_discard_reads option"
"$meta_executable" \
    --paired \
    --input "$test_dir/scrb_seq_fastq.1.gz" \
    --read2_in "$test_dir/scrb_seq_fastq.2.gz" \
    --bc_pattern CCCCCCNNNNNNNNNN \
    --bc_pattern2 CCCCCCNNNNNNNNNN \
    --umitools_extract_method string \
    --umitools_umi_separator '_' \
    --umitools_grouping_method directional \
    --umi_discard_read 2 \
    --read1_out scrb_seq_fastq.1.umi_extract.fastq.gz \
    --random_seed 1

echo ">> Checking if the correct files are present"
[ ! -f "scrb_seq_fastq.1.umi_extract.fastq.gz" ] && echo "Read 1 file is missing" && exit 1
[ ! -s "scrb_seq_fastq.1.umi_extract.fastq.gz" ] && echo "Read 1 file is empty" && exit 1
[ -f "scrb_seq_fastq.2.umi_extract.fastq.gz" ] && echo "Read 2 is not discarded" && exit 1

echo ">> Checking if the files are correct"
diff -q <(gunzip -c "${meta_resources_dir}/scrb_seq_fastq.1.umi_extract.fastq.gz") <(gunzip -c "$test_dir/scrb_seq_fastq.1.umi_extract.fastq.gz") || \
    (echo "Read 1 file is not correct" && exit 1)

rm scrb_seq_fastq.1.umi_extract.fastq.gz

############################################################################################################

echo ">>> Test 3: Testing for single-end reads"
"$meta_executable" \
    --input "$test_dir/slim.fastq.gz" \
    --bc_pattern "^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})" \
    --umitools_extract_method regex \
    --umitools_umi_separator '_' \
    --umitools_grouping_method directional \
    --read1_out slim.umi_extract.fastq.gz \
    --random_seed 1 

echo ">> Checking if the correct files are present"
[ ! -f "slim.umi_extract.fastq.gz" ] && echo "Trimmed reads file missing" && exit 1
[ ! -s "slim.umi_extract.fastq.gz" ] && echo "Trimmed reads file is empty" && exit 1

echo ">> Checking if the files are correct"
diff -q <(gunzip -c "${meta_resources_dir}/slim.umi_extract.fastq.gz") <(gunzip -c "$test_dir/slim.umi_extract.fastq.gz") || \
    (echo "Trimmed reads file is not correct" && exit 1)

rm slim.umi_extract.fastq.gz

echo ">>> Test finished successfully"
exit 0