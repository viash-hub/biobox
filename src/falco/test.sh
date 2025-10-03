#!/bin/bash

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Prepare test data"
echo ">> Preparing test data"

test_data_dir="$meta_temp_dir/test_data"
mkdir $test_data_dir

create_test_fastq $test_data_dir/R1.fastq.gz 200 50
create_test_fastq $test_data_dir/R2.fastq.gz 200 50

log "Run falco on test data, output to dir"
test_1_output="$meta_temp_dir/output1"
$meta_executable \
  --input "$test_data_dir/R1.fastq.gz;$test_data_dir/R2.fastq.gz" \
  --outdir "$test_1_output"

log "Checking whether output exists"
check_dir_exists "$test_1_output" "output directory (when outputting to directory)"
check_file_exists "$test_1_output/R1.fastq.gz_fastqc_report.html" "report for R1 (when outputting to directory)"
check_file_not_empty "$test_1_output/R1.fastq.gz_fastqc_report.html" "report for R1 (when outputting to directory)"
check_file_exists "$test_1_output/R1.fastq.gz_summary.txt" "summary for R1 (when outputting to directory)"
check_file_not_empty "$test_1_output/R1.fastq.gz_summary.txt" "summary for R1 (when outputting to directory)"
check_file_exists "$test_1_output/R1.fastq.gz_fastqc_data.txt" "data for R1 (when outputting to directory)"
check_file_not_empty "$test_1_output/R1.fastq.gz_fastqc_data.txt" "data for R1 (when outputting to directory)"
check_file_exists "$test_1_output/R2.fastq.gz_fastqc_report.html" "report for R2 (when outputting to directory)"
check_file_not_empty "$test_1_output/R2.fastq.gz_fastqc_report.html" "report for R2 (when outputting to directory)"
check_file_exists "$test_1_output/R2.fastq.gz_summary.txt" "summary for R2 (when outputting to directory)"
check_file_not_empty "$test_1_output/R2.fastq.gz_summary.txt" "summary for R2 (when outputting to directory)"
check_file_exists "$test_1_output/R2.fastq.gz_fastqc_data.txt" "data for R2 (when outputting to directory)"
check_file_not_empty "$test_1_output/R2.fastq.gz_fastqc_data.txt" "data for R2 (when outputting to directory)"

log "cleanup"
rm -rf "$test_1_output"

log "Run falco on test data, output to individual files. Please note this is only possible for 1 input fastq file!"
test_2_output="$meta_temp_dir/output2/"
$meta_executable \
  --input "$test_data_dir/R1.fastq.gz" \
  --data_filename "$test_2_output/data.txt" \
  --report_filename "$test_2_output/report.html" \
  --summary_filename "$test_2_output/summary.txt" \
  --outdir "$test_2_output"

log "Checking whether output exists"
check_dir_exists "$test_2_output" "subsampling test output"
check_file_exists "$test_2_output/report.html" "report output (after subsampling)"
check_file_not_empty "$test_2_output/report.html" "report output (after subsampling)"
check_file_exists "$test_2_output/summary.txt" "summary output (after subsampling)"
check_file_not_empty "$test_2_output/summary.txt" "summary output (after subsampling)"
check_file_exists "$test_2_output/data.txt" "data output (after subsampling)"
check_file_not_empty "$test_2_output/data.txt" "data output (after subsampling)"

log "Cleanup"
rm -rf "$test_2_output"

echo "Run falco on test data, subsample"
test_3_output="$meta_temp_dir/output3/"

$meta_executable \
  --input "$test_data_dir/R1.fastq.gz" \
  --data_filename "$test_3_output/data.txt" \
  --report_filename "$test_3_output/report.html" \
  --summary_filename "$test_3_output/summary.txt" \
  --subsample 100 \
  --outdir "$test_3_output"

log "Checking whether output exists"
check_dir_exists "$test_3_output" "subsampling test output"
check_file_exists "$test_3_output/report.html" "report output (after subsampling)"
check_file_not_empty "$test_3_output/report.html" "report output (after subsampling)"
check_file_exists "$test_3_output/summary.txt" "summary output (after subsampling)"
check_file_not_empty "$test_3_output/summary.txt" "summary output (after subsampling)"
check_file_exists "$test_3_output/data.txt" "data output (after subsampling)"
check_file_not_empty "$test_3_output/data.txt" "data output (after subsampling)"

log "cleanup"
rm -rf "$meta_temp_dir/output3/"

log "Run falco with empty fastq"
touch $test_data_dir/empty_R1.fastq
test_4_output="$meta_temp_dir/output4/"

$meta_executable \
  --input "$test_data_dir/empty_R1.fastq" \
  --data_filename "$test_4_output/data.txt" \
  --report_filename "$test_4_output/report.html" \
  --summary_filename "$test_4_output/summary.txt" \
  --allow_empty_input \
  --outdir "$test_4_output"

log "Checking whether output exists"

check_dir_exists "$test_4_output" "output directory (empty input files)"
check_file_exists "$test_4_output/report.html" "report output (empty input files)"
check_file_not_empty "$test_4_output/report.html" "report output (empty input files)"
check_file_exists "$test_4_output/summary.txt" "summary output (empty input files)"
check_file_not_empty "$test_4_output/summary.txt" "summary output (empty input files)"
check_file_exists "$test_4_output/data.txt" "data output (empty input files)"
check_file_not_empty "$test_4_output/data.txt" "data output (empty input files)"

log "cleanup"
rm -rf "$test_4_output"

print_test_summary