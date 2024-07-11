#!/bin/bash

# Download test data
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/slim.fastq.gz
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/scrb_seq_fastq.1.gz
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/scrb_seq_fastq.2.gz

gunzip -f slim.fastq.gz scrb_seq_fastq.1.gz scrb_seq_fastq.2.gz

# smaller datasets
head -n 120 slim.fastq > slim_30.fastq
head -n 120 scrb_seq_fastq.1 > scrb_seq_fastq.1_30
head -n 120 scrb_seq_fastq.2 > scrb_seq_fastq.2_30
rm slim.fastq scrb_seq_fastq.1 scrb_seq_fastq.2

# Generate expected output
# Test 1 and 2
umi_tools extract \
    --stdin "scrb_seq_fastq.1_30" \
    --read2-in "scrb_seq_fastq.2_30" \
    --bc-pattern "CCCCCCNNNNNNNNNN" \
    --bc-pattern2 "CCCCCCNNNNNNNNNN" \
    --extract-method string \
    --stdout scrb_seq_fastq.1_30.extract \
    --read2-out scrb_seq_fastq.2_30.extract \
    --random-seed 1

# Test 3
umi_tools extract \
    --stdin "slim_30.fastq" \
    --bc-pattern "^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})" \
    --extract-method regex \
    --stdout slim_30.extract \
    --random-seed 1