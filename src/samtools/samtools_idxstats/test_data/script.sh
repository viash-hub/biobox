#!/bin/bash

# dowload test data from snakemake wrapper
if [ ! -d /tmp/idxstats_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers.git /tmp/idxstats_source
fi

cp -r /tmp/idxstats_source/bio/samtools/idxstats/test/mapped/* src/samtools/idxstats/test_data
# samtools idxstats a.sorted.bam > a.sorted.idxstats

# dowload test data from nf-core module
wget https://github.com/nf-core/test-datasets/raw/modules/data/genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam
wget https://github.com/nf-core/test-datasets/raw/modules/data/genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam.bai
# samtools idxstats test.paired_end.sorted.bam > test_ref.paired_end.sorted.idxstats