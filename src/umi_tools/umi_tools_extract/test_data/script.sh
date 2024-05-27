#!/bin/bash

# Download test data
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/slim.fastq.gz
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/scrb_seq_fastq.1.gz
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/scrb_seq_fastq.2.gz

# unzip
gunzip scrb_seq_fastq.1.gz
gunzip scrb_seq_fastq.2.gz
gunzip slim.fastq.gz

# get a subsample of the reads (50 first reads from each file)
head -n 100 scrb_seq_fastq.1 > scrb_seq_fastq.1.subsample.fastq
head -n 100 scrb_seq_fastq.2 > scrb_seq_fastq.2.subsample.fastq
head -n 200 slim.fastq > slim.subsample.fastq

rm scrb_seq_fastq.1 scrb_seq_fastq.2 slim.fastq
