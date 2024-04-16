#!/bin/bash

# Tests are sourced from:
#   https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-direct-demultiplexing-bcl-convert
# Test input files are fetched from:
# https://cf.10xgenomics.com/supp/spatial-exp/demultiplexing/iseq-DI.tar.gz
# https://cf.10xgenomics.com/supp/spatial-exp/demultiplexing/bcl_convert_samplesheet.csv

set -eo pipefail

echo ">> Fetching and preparing test data"
data_src="https://cf.10xgenomics.com/supp/spatial-exp/demultiplexing/iseq-DI.tar.gz"
sample_sheet_src="https://cf.10xgenomics.com/supp/spatial-exp/demultiplexing/bcl_convert_samplesheet.csv"
test_data_dir="test_data"

mkdir $test_data_dir
wget -q $data_src -O $test_data_dir/data.tar.gz
wget -q $sample_sheet_src -O $test_data_dir/sample_sheet.csv
tar xzf $test_data_dir/data.tar.gz -C $test_data_dir
rm $test_data_dir/data.tar.gz

echo ">> Execute and verify output"

$meta_executable \
  --bcl_input_directory "$test_data_dir/iseq-DI" \
  --sample_sheet "$test_data_dir/sample_sheet.csv" \
  --output_directory fastq \
  --reports reports \
  --logs logs

echo ">>> Checking whether the output dir exists"
[[ ! -d fastq ]] && echo "Output dir could not be found!" && exit 1

echo ">>> Checking whether output fastq files are created"
[[ ! -f fastq/Undetermined_S0_L001_R1_001.fastq.gz ]] && echo "Output fastq files could not be found!" && exit 1
[[ ! -f fastq/iseq-DI_S1_L001_R1_001.fastq.gz ]] && echo "Output fastq files could not be found!" && exit 1

echo ">>> Checking whether the report dir exists"
[[ ! -d reports ]] && echo "Reports dir could not be found!" && exit 1

echo ">>> Checking whether the log dir exists"
[[ ! -d logs ]] && echo "Logs dir could not be found!" && exit 1

# print final message
echo ">>> Test finished successfully"

echo ">> Execute with additional arguments and verify output"

$meta_executable \
  --bcl_input_directory "$test_data_dir/iseq-DI" \
  --sample_sheet "$test_data_dir/sample_sheet.csv" \
  --output_directory fastq1 \
  --bcl_only_matched_reads true \
  --bcl_num_compression_threads 1 \
  --no_lane_splitting false \
  --fastq_gzip_compression_level 9

echo ">> Checking whether the output dir exists"
[[ ! -d fastq1 ]] && echo "Output dir could not be found!" && exit 1

echo ">> Checking whether output fastq files are created"
[[ -f fastq1/Undetermined_S0_L001_R1_001.fastq.gz ]] && echo "Undetermined should not be generated!" && exit 1
[[ ! -f fastq1/iseq-DI_S1_L001_R1_001.fastq.gz ]] && echo "Output fastq files could not be found!" && exit 1

# print final message
echo ">> Test finished successfully"

# do not remove this
# as otherwise your test might exit with a different exit code
exit 0
