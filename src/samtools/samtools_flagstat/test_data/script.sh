#!/bin/bash

# Download test data from snakemake wrapper

wget https://raw.githubusercontent.com/snakemake/snakemake-wrappers/3a4f7004281efc176fd9af732ad88d00c47d432d/bio/samtools/flagstat/test/mapped/a.bam
samtools index a.bam
# samtools flagstat a.bam > a_ref.flagstat


# Download test data from nf-core module

wget https://github.com/nf-core/test-datasets/raw/modules/data/genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam
wget https://github.com/nf-core/test-datasets/raw/modules/data/genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam.bai
# samtools flagstat test.paired_end.sorted.bam > test_ref.paired_end.sorted.flagstat