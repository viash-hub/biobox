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

log "Generating test reference genome..."
create_test_fasta "$test_data_dir/test_ref.fasta" 2 1000
check_file_exists "$test_data_dir/test_ref.fasta" "test reference genome"

# --- Test Case 1: Basic indexing ---
log "Starting TEST 1: Basic minibwa indexing"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/minibwa_index"

log "Validating TEST 1 outputs..."
check_dir_exists "$meta_temp_dir/minibwa_index" "output index directory"

index_files=(
  "$meta_temp_dir/minibwa_index/test_ref.l2b"
  "$meta_temp_dir/minibwa_index/test_ref.mbw"
)

for file in "${index_files[@]}"; do
  check_file_exists "$file" "minibwa index file $(basename "$file")"
  check_file_not_empty "$file" "minibwa index file $(basename "$file")"
done

check_file_not_exists "$meta_temp_dir/minibwa_index/test_ref.meth.mbw" "bisulfite index file (should not be built by default)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Custom prefix ---
log "Starting TEST 2: minibwa indexing with custom prefix"

log "Executing $meta_name with custom prefix..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/custom_index" \
  --prefix "custom_genome"

log "Validating TEST 2 outputs..."
check_dir_exists "$meta_temp_dir/custom_index" "custom index directory"

custom_index_files=(
  "$meta_temp_dir/custom_index/custom_genome.l2b"
  "$meta_temp_dir/custom_index/custom_genome.mbw"
)

for file in "${custom_index_files[@]}"; do
  check_file_exists "$file" "custom-prefixed index file $(basename "$file")"
  check_file_not_empty "$file" "custom-prefixed index file $(basename "$file")"
done

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Bisulfite index ---
log "Starting TEST 3: minibwa indexing with bisulfite (--meth) index"

log "Executing $meta_name with --meth..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/bisulfite_index" \
  --meth

log "Validating TEST 3 outputs..."
check_dir_exists "$meta_temp_dir/bisulfite_index" "bisulfite index directory"

bisulfite_index_files=(
  "$meta_temp_dir/bisulfite_index/test_ref.l2b"
  "$meta_temp_dir/bisulfite_index/test_ref.mbw"
  "$meta_temp_dir/bisulfite_index/test_ref.meth.mbw"
)

for file in "${bisulfite_index_files[@]}"; do
  check_file_exists "$file" "bisulfite index file $(basename "$file")"
  check_file_not_empty "$file" "bisulfite index file $(basename "$file")"
done

log "✅ TEST 3 completed successfully"

print_test_summary "All tests"
