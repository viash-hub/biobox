#!/bin/bash

# dowload test data from snakemake wrapper
if [ ! -d /tmp/fastq_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers.git /tmp/fastq_source
fi

cp -r /tmp/fastq_source/bio/samtools/fastx/test/*.sam src/samtools/samtools_fastq/test_data/