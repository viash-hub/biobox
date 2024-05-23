#!/bin/bash

# dowload test data from snakemake wrapper
if [ ! -d /tmp/view_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers.git /tmp/view_source
fi

cp -r /tmp/idxstats_source/bio/samtools/view/test/*.sam src/samtools/samtools_view/test_data