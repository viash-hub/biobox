#!/bin/bash

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Starting tests for $meta_name"

test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

create_test_fastq "$test_data_dir/read1_rep1.fastq" 4 35
create_test_fastq "$test_data_dir/read1_rep2.fastq" 6 35
create_test_fastq "$test_data_dir/read2_rep1.fastq" 4 35
create_test_fastq "$test_data_dir/read2_rep2.fastq" 6 35

gzip -c "$test_data_dir/read1_rep1.fastq" > "$test_data_dir/read1_rep1.fastq.gz"
gzip -c "$test_data_dir/read1_rep2.fastq" > "$test_data_dir/read1_rep2.fastq.gz"
gzip -c "$test_data_dir/read2_rep1.fastq" > "$test_data_dir/read2_rep1.fastq.gz"
gzip -c "$test_data_dir/read2_rep2.fastq" > "$test_data_dir/read2_rep2.fastq.gz"

log "Starting TEST 1: paired-end, multiple replicates, uncompressed input"
out_dir="$meta_temp_dir/output1"
mkdir -p "$out_dir"
"$meta_executable" \
  --read_1 "$test_data_dir/read1_rep1.fastq;$test_data_dir/read1_rep2.fastq" \
  --read_2 "$test_data_dir/read2_rep1.fastq;$test_data_dir/read2_rep2.fastq" \
  --fastq_1 "$out_dir/read_1.merged.fastq" \
  --fastq_2 "$out_dir/read_2.merged.fastq"

check_file_exists "$out_dir/read_1.merged.fastq" "merged read 1 (uncompressed)"
check_file_not_empty "$out_dir/read_1.merged.fastq" "merged read 1 (uncompressed)"
check_file_line_count "$out_dir/read_1.merged.fastq" 40 "merged read 1 (uncompressed)"
check_file_exists "$out_dir/read_2.merged.fastq" "merged read 2 (uncompressed)"
check_file_not_empty "$out_dir/read_2.merged.fastq" "merged read 2 (uncompressed)"
check_file_line_count "$out_dir/read_2.merged.fastq" 40 "merged read 2 (uncompressed)"
log "TEST 1 completed successfully"

log "Starting TEST 2: paired-end, multiple replicates, gzipped input"
out_dir="$meta_temp_dir/output2"
mkdir -p "$out_dir"
"$meta_executable" \
  --read_1 "$test_data_dir/read1_rep1.fastq.gz;$test_data_dir/read1_rep2.fastq.gz" \
  --read_2 "$test_data_dir/read2_rep1.fastq.gz;$test_data_dir/read2_rep2.fastq.gz" \
  --fastq_1 "$out_dir/read_1.merged.fastq" \
  --fastq_2 "$out_dir/read_2.merged.fastq"

check_file_exists "$out_dir/read_1.merged.fastq" "merged read 1 (gzipped input)"
check_file_line_count "$out_dir/read_1.merged.fastq" 40 "merged read 1 (gzipped input)"
check_file_exists "$out_dir/read_2.merged.fastq" "merged read 2 (gzipped input)"
check_file_line_count "$out_dir/read_2.merged.fastq" 40 "merged read 2 (gzipped input)"
log "TEST 2 completed successfully"

log "Starting TEST 3: paired-end, mixed compressed/uncompressed replicates"
out_dir="$meta_temp_dir/output3_mixed"
mkdir -p "$out_dir"
"$meta_executable" \
  --read_1 "$test_data_dir/read1_rep1.fastq;$test_data_dir/read1_rep2.fastq.gz" \
  --read_2 "$test_data_dir/read2_rep1.fastq.gz;$test_data_dir/read2_rep2.fastq" \
  --fastq_1 "$out_dir/read_1.merged.fastq" \
  --fastq_2 "$out_dir/read_2.merged.fastq"

check_file_exists "$out_dir/read_1.merged.fastq" "merged read 1 (mixed compression)"
check_file_line_count "$out_dir/read_1.merged.fastq" 40 "merged read 1 (mixed compression)"
check_file_exists "$out_dir/read_2.merged.fastq" "merged read 2 (mixed compression)"
check_file_line_count "$out_dir/read_2.merged.fastq" 40 "merged read 2 (mixed compression)"
log "TEST 3 completed successfully"

log "Starting TEST 4: single-end, multiple replicates"
out_dir="$meta_temp_dir/output4"
mkdir -p "$out_dir"
"$meta_executable" \
  --read_1 "$test_data_dir/read1_rep1.fastq;$test_data_dir/read1_rep2.fastq" \
  --fastq_1 "$out_dir/read_1.merged.fastq"

check_file_exists "$out_dir/read_1.merged.fastq" "merged read 1 (single-end)"
check_file_line_count "$out_dir/read_1.merged.fastq" 40 "merged read 1 (single-end)"
check_file_not_exists "$out_dir/read_2.merged.fastq" "read 2 output (single-end, should not be created)"
log "TEST 4 completed successfully"

print_test_summary "cat_fastq"
