#!/bin/bash

# dowload test data from snakemake wrapper
if [ ! -d /tmp/idxstats_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers.git /tmp/idxstats_source
fi

cp -r /tmp/idxstats_source/bio/samtools/idxstats/test/mapped/* src/samtools/idxstats/test_data
# samtools index a_ref.sorted.bam -o a_ref.sorted.bam.bai
# samtools index a_ref.sorted.bam -c a_ref.sorted.bam.csi


