#!/bin/bash

# dowload test data from snakemake wrapper
if [ ! -d /tmp/idxstats_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers.git /tmp/sort_source
fi

cp -r /tmp/sort_source/bio/samtools/sort/test/mapped/* src/samtools/samtools_sort/test_data
