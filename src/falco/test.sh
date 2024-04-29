#!/bin/bash

set -e

echo "> Prepare test data"

# We use data from this repo: https://github.com/hartwigmedical/testData
echo ">> Fetching and preparing test data"
fastq1="https://github.com/hartwigmedical/testdata/raw/master/100k_reads_hiseq/TESTX/TESTX_H7YRLADXX_S1_L001_R1_001.fastq.gz"
fastq2="https://github.com/hartwigmedical/testdata/raw/master/100k_reads_hiseq/TESTX/TESTX_H7YRLADXX_S1_L001_R2_001.fastq.gz"
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

test_data_dir="$TMPDIR/test_data"

mkdir $test_data_dir
wget -q $fastq1 -O $test_data_dir/R1.fastq.gz
wget -q $fastq2 -O $test_data_dir/R2.fastq.gz

echo ">> Run falco on test data, output to dir"
echo ">>> Run falco"
$meta_executable \
  --input "$test_data_dir/R1.fastq.gz;$test_data_dir/R2.fastq.gz" \
  --outdir out/

echo ">>> Checking whether output exists"
[ ! -d "out" ] && echo "Output directory not created" && exit 1
[ ! -f "out/R1.fastq.gz_fastqc_report.html" ] && echo "Report not created" && exit 1
[ ! -f "out/R1.fastq.gz_summary.txt" ] && echo "Summary not created" && exit 1
[ ! -f "out/R1.fastq.gz_fastqc_data.txt" ] && echo "fastqc_data not created" && exit 1
[ ! -f "out/R2.fastq.gz_fastqc_report.html" ] && echo "Report not created" && exit 1
[ ! -f "out/R2.fastq.gz_summary.txt" ] && echo "Summary not created" && exit 1
[ ! -f "out/R2.fastq.gz_fastqc_data.txt" ] && echo "fastqc_data not created" && exit 1

echo ">>> cleanup"
rm -rf out/

echo ">> Run falco on test data, output to individual files"
echo ">>> Please note this is only possible for 1 input fastq file!"
echo ">>> Run falco"
$meta_executable \
  --input "$test_data_dir/R1.fastq.gz" \
  --data_filename out/data.txt \
  --report_filename out/report.html \
  --summary_filename out/summary.txt \
  --outdir out/

echo ">>> Checking whether output exists"
[ ! -d "out" ] && echo "Output directory not created" && exit 1
[ ! -f "out/report.html" ] && echo "Report not created" && exit 1
[ ! -f "out/summary.txt" ] && echo "Summary not created" && exit 1
[ ! -f "out/data.txt" ] && echo "fastqc_data not created" && exit 1

echo ">>> cleanup"
rm -rf out/

echo ">> Run falco on test data, subsample"
echo ">>> Run falco"
$meta_executable \
  --input "$test_data_dir/R1.fastq.gz" \
  --data_filename out/data.txt \
  --report_filename out/report.html \
  --summary_filename out/summary.txt \
  --subsample 100 \
  --outdir out

echo ">>> Checking whether output exists"
[ ! -d "out" ] && echo "Output directory not created" && exit 1
[ ! -f "out/report.html" ] && echo "Report not created" && exit 1
[ ! -f "out/summary.txt" ] && echo "Summary not created" && exit 1
[ ! -f "out/data.txt" ] && echo "fastqc_data not created" && exit 1

echo ">>> cleanup"
rm -rf out/

echo "All tests succeeded!"
