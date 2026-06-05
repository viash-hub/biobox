#!/bin/bash

## VIASH START
## VIASH END

# source the helpers & setup test env
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Test for $meta_name starting ..."

# create and validate test-data
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir" "$meta_temp_dir/index"
create_test_fasta "$test_data_dir/22_20-21M.fa" 1 10000
check_file_exists "$test_data_dir/22_20-21M.fa" "test reference FASTA"

# --- TEST: create a hisat2 index ---
log "Create a hisat2 index from a reference FASTA"
"$meta_executable" \
  --reference "$test_data_dir/22_20-21M.fa" \
  --index_dir "$meta_temp_dir/index" \
  --index_prefix genome

for ext in 1.ht2 2.ht2 3.ht2 4.ht2 5.ht2 6.ht2 7.ht2 8.ht2; do
  check_file_exists "$meta_temp_dir/index/genome.$ext" "index file genome.$ext"
done

log "✓ TEST completed successfully: 8 index files were created"

print_test_summary "$meta_name tests"