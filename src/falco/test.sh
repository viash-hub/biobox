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
  --outdir "$TMPDIR/output1"

echo ">>> Checking whether output exists"
[ ! -d "$TMPDIR/output1" ] && echo "Output directory not created" && exit 1
[ ! -f "$TMPDIR/output1/R1.fastq.gz_fastqc_report.html" ] && echo "Report not created" && exit 1
[ ! -f "$TMPDIR/output1/R1.fastq.gz_summary.txt" ] && echo "Summary not created" && exit 1
[ ! -f "$TMPDIR/output1/R1.fastq.gz_fastqc_data.txt" ] && echo "fastqc_data not created" && exit 1
[ ! -f "$TMPDIR/output1/R2.fastq.gz_fastqc_report.html" ] && echo "Report not created" && exit 1
[ ! -f "$TMPDIR/output1/R2.fastq.gz_summary.txt" ] && echo "Summary not created" && exit 1
[ ! -f "$TMPDIR/output1/R2.fastq.gz_fastqc_data.txt" ] && echo "fastqc_data not created" && exit 1

echo ">>> cleanup"
rm -rf "$TMPDIR/output1"

echo ">> Run falco on test data, output to individual files"
echo ">>> Please note this is only possible for 1 input fastq file!"
echo ">>> Run falco"
$meta_executable \
  --input "$test_data_dir/R1.fastq.gz" \
  --data_filename "$TMPDIR/output2/data.txt" \
  --report_filename "$TMPDIR/output2/report.html" \
  --summary_filename "$TMPDIR/output2/summary.txt" \
  --outdir "$TMPDIR/output2/"

echo ">>> Checking whether output exists"
[ ! -d "$TMPDIR/output2" ] && echo "Output directory not created" && exit 1
[ ! -f "$TMPDIR/output2/report.html" ] && echo "Report not created" && exit 1
[ ! -f "$TMPDIR/output2/summary.txt" ] && echo "Summary not created" && exit 1
[ ! -f "$TMPDIR/output2/data.txt" ] && echo "fastqc_data not created" && exit 1

echo ">>> cleanup"
rm -rf $TMPDIR/output2/

echo ">> Run falco on test data, subsample"
echo ">>> Run falco"
$meta_executable \
  --input "$test_data_dir/R1.fastq.gz" \
  --data_filename "$TMPDIR/output3/data.txt" \
  --report_filename "$TMPDIR/output3/report.html" \
  --summary_filename "$TMPDIR/output3/summary.txt" \
  --subsample 100 \
  --outdir "$TMPDIR/output3"

echo ">>> Checking whether output exists"
[ ! -d "$TMPDIR/output3" ] && echo "Output directory not created" && exit 1
[ ! -f "$TMPDIR/output3/report.html" ] && echo "Report not created" && exit 1
[ ! -f "$TMPDIR/output3/summary.txt" ] && echo "Summary not created" && exit 1
[ ! -f "$TMPDIR/output3/data.txt" ] && echo "fastqc_data not created" && exit 1

echo ">>> cleanup"
rm -rf "$TMPDIR/output3/"

echo "All tests succeeded!"
