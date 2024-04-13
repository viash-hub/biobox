#!/bin/bash

# dowload test data from nf-core module
wget https://github.com/nf-core/test-datasets/raw/modules/data/genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam
wget https://github.com/nf-core/test-datasets/raw/modules/data/genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam.bai
# samtools stats test.paired_end.sorted.bam > ref.paired_end.sorted.txt