#!/bin/bash

## VIASH START
## VIASH END

# source the helpers
source "$meta_resources_dir/test_helpers.sh"
# setup test env
setup_test_env

log "Starting tests for $meta_name"

# # donwloading test sample data from hisat2 repo
# curl -L https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.2.tar.gz \
#     | tar -xz -C "$test_dir" \
#     --strip-components=3 \
#     hisat2-2.2.2/example/reference/22_20-21M.fa

# create and validate test-data
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"
create_test_fasta "$test_data_dir/22_20-21M.fa" 1 10000
check_file_exists "$test_data_dir/22_20-21M.fa" "test reference FASTA"

# --- TEST: create a hisat2 index ---
log "Create a hisat2 index from a reference FASTA"
"$meta_executable" \
  --input "$test_data_dir/22_20-21M.fa" \
  --index_dir "$meta_temp_dir/index" \
  --index_prefix genome

log "✓ TEST completed successfully: 8 index files were created"

print_test_summary "$meta_name tests"