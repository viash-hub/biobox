#!/bin/bash

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/rsem.tar.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357070_1.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357070_2.fastq.gz

# decompress file without keeping the original
gunzip SRR6357070_1.fastq.gz
gunzip SRR6357070_2.fastq.gz

# only keep 100 reads per file
head -n 400 SRR6357070_1.fastq > SRR6357070_1.fastq.tmp
head -n 400 SRR6357070_2.fastq > SRR6357070_2.fastq.tmp

mv SRR6357070_1.fastq.tmp SRR6357070_1.fastq
mv SRR6357070_2.fastq.tmp SRR6357070_2.fastq

rm SRR6357070_1.fastq.gz
rm SRR6357070_2.fastq.gz

gzip SRR6357070_1.fastq
gzip SRR6357070_2.fastq

