#!/bin/bash

## VIASH START
## VIASH END

# source the helpers
source "$meta_resources_dir/test_helpers.sh"
# setup test env
setup_test_env

log "Starting tests for $meta_name"

# --- TEST_1: basic functionality ---
log "Starting TEST_1: basic functionality"

# --- create and validate test-data ---
log "Creating & fetching test data"

test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir/test_1"

# create the fastq read
create_test_fastq "$test_data_dir/test_1/read_1.fastq" 1 10000
check_file_exists "$test_data_dir/test_1/read_1.fastq" "test FASTQ file"

# opted here to fetch data from the HISAT2 repo because of the specific file format (fasta reads)
# prefer it to not be dependent on the hisat2-build command for isolation
curl -kL https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.2.tar.gz \
    | tar -xz -C "$test_data_dir/test_1/" \
    --strip-components=2 \
    hisat2-2.2.2/example/index \
    hisat2-2.2.2/example/reads   
# check if fetched data exists
for i in $(seq 1 8); do
  check_file_exists "$test_data_dir/test_1/index/22_20-21M_snp.${i}.ht2" "HISAT2 index file ${i}"
done

check_file_exists "$test_data_dir/test_1/reads/reads_1.fa" "reads file 1"
check_file_exists "$test_data_dir/test_1/reads/reads_2.fa" "reads file 2"

# --- actual testing ---
log "Test data created & fetched, now executing $meta_name with basic parameters..."
"$meta_executable" \
  --index_dir "$test_data_dir/test_1/index" \
  --index_prefix 22_20-21M_snp \
  --input "$test_data_dir/test_1/read_1.fastq" \
  --output_sam "$test_data_dir/test1/output.sam"

check_file_exists "$test_data_dir/test1/output.sam" "Output SAM file"
check_file_not_empty "$test_data_dir/test1/output.sam" "Output SAM file"

log "✓ TEST_1 completed successfully"

#TODO: test for fasta flag

# (# test for alignment scores)
