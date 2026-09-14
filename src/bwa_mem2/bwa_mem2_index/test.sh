#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# Create test data directory
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

# --- Test Case 1: Basic indexing ---
log "Starting TEST 1: Basic BWA indexing"

log "Generating test reference genome..."
create_test_fasta "$test_data_dir/test_ref.fasta" 2 1000
check_file_exists "$test_data_dir/test_ref.fasta" "test reference genome"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/bwa_index"

log "Validating TEST 1 outputs..."
check_dir_exists "$meta_temp_dir/bwa_index" "output index directory"

# Check for BWA index files with the prefix used by bwa index
index_files=(
  "$meta_temp_dir/bwa_index/test_ref.0123"
  "$meta_temp_dir/bwa_index/test_ref.amb"
  "$meta_temp_dir/bwa_index/test_ref.ann"
  "$meta_temp_dir/bwa_index/test_ref.bwt.2bit.64"
  "$meta_temp_dir/bwa_index/test_ref.pac"
)

for file in "${index_files[@]}"; do
  check_file_exists "$file" "BWA index file $(basename "$file")"
  check_file_not_empty "$file" "BWA index file $(basename "$file")"
done

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Custom prefix ---
log "Starting TEST 2: BWA indexing with custom prefix"

log "Executing $meta_name with custom prefix..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/custom_index" \
  --prefix "custom_genome"

log "Validating TEST 2 outputs..."
check_dir_exists "$meta_temp_dir/custom_index" "custom index directory"

# Check for index files with custom prefix
log "Checking for custom-prefixed index files..."
custom_index_files=(
  "$meta_temp_dir/custom_index/custom_genome.0123"
  "$meta_temp_dir/custom_index/custom_genome.amb"
  "$meta_temp_dir/custom_index/custom_genome.ann"
  "$meta_temp_dir/custom_index/custom_genome.bwt.2bit.64"
  "$meta_temp_dir/custom_index/custom_genome.pac"
)

for file in "${custom_index_files[@]}"; do
  check_file_exists "$file" "custom-prefixed index file $(basename "$file")"
  check_file_not_empty "$file" "custom-prefixed index file $(basename "$file")"
done

log "✅ TEST 2 completed successfully"

print_test_summary "All tests completed successfully"
