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

# --- Test Case 1: Default indexing ---
log "Starting TEST 1: Default BWA indexing"

log "Generating test reference genome..."
create_test_fasta "$test_data_dir/test_ref.fasta" 2 1000
check_file_exists "$test_data_dir/test_ref.fasta" "test reference genome"

log "Executing $meta_name with default parameters..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/bwa_index"

log "Validating TEST 1 outputs..."
check_dir_exists "$meta_temp_dir/bwa_index" "output index directory"

# Check for index files with the default prefix
index_files=(
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

# --- Test Case 3: Methylation parameter ---
log "Starting TEST 3: BWA indexing with methylation parameter"

log "Executing $meta_name with methylation parameter..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/methylation_index" \
  --meth

log "Validating TEST 3 outputs..."
check_dir_exists "$meta_temp_dir/methylation_index" "methylation index directory"

# Check for normal and methylation seed index files (default prefix)
log "Checking for normal and methylation seed index files..."
meth_index_files=(
  "$meta_temp_dir/methylation_index/test_ref.amb"
  "$meta_temp_dir/methylation_index/test_ref.ann"
  "$meta_temp_dir/methylation_index/test_ref.bwt.2bit.64"
  "$meta_temp_dir/methylation_index/test_ref.pac"
  "$meta_temp_dir/methylation_index/test_ref.meth.amb"
  "$meta_temp_dir/methylation_index/test_ref.meth.ann"
  "$meta_temp_dir/methylation_index/test_ref.meth.bwt.2bit.64"
  "$meta_temp_dir/methylation_index/test_ref.meth.pac"
  "$meta_temp_dir/methylation_index/test_ref.meth.fa"
)

for file in "${meth_index_files[@]}"; do
  check_file_exists "$file" "methylation index file $(basename "$file")"
  check_file_not_empty "$file" "methylation index file $(basename "$file")"
done

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Unpacked reference ---
log "Starting TEST 4: BWA indexing with unpacked reference"

log "Executing $meta_name with --emit_unpacked_ref..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/unpacked_ref_index" \
  --emit_unpacked_ref

log "Validating TEST 4 outputs..."
check_dir_exists "$meta_temp_dir/unpacked_ref_index" "unpacked reference index directory"

# Check for index files plus the unpacked .0123 reference (default prefix)
log "Checking for index files and unpacked reference..."
unpacked_index_files=(
  "$meta_temp_dir/unpacked_ref_index/test_ref.0123"
  "$meta_temp_dir/unpacked_ref_index/test_ref.amb"
  "$meta_temp_dir/unpacked_ref_index/test_ref.ann"
  "$meta_temp_dir/unpacked_ref_index/test_ref.bwt.2bit.64"
  "$meta_temp_dir/unpacked_ref_index/test_ref.pac"
)

for file in "${unpacked_index_files[@]}"; do
  check_file_exists "$file" "unpacked reference index file $(basename "$file")"
  check_file_not_empty "$file" "unpacked reference index file $(basename "$file")"
done

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Suffix-array sample rate ---
log "Starting TEST 5: BWA indexing with custom suffix-array sample rate"

log "Executing $meta_name with --sa_sample_rate 1..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/sa_rate_index" \
  --sa_sample_rate 1

log "Validating TEST 5 outputs..."
check_dir_exists "$meta_temp_dir/sa_rate_index" "sample rate index directory"

sa_rate_index_files=(
  "$meta_temp_dir/sa_rate_index/test_ref.amb"
  "$meta_temp_dir/sa_rate_index/test_ref.ann"
  "$meta_temp_dir/sa_rate_index/test_ref.bwt.2bit.64"
  "$meta_temp_dir/sa_rate_index/test_ref.pac"
)

for file in "${sa_rate_index_files[@]}"; do
  check_file_exists "$file" "sample rate index file $(basename "$file")"
  check_file_not_empty "$file" "sample rate index file $(basename "$file")"
done

# A denser suffix array (-u 1) must produce a larger index than the default (-u 3, from TEST 1)
default_size=$(wc -c < "$meta_temp_dir/bwa_index/test_ref.bwt.2bit.64")
dense_size=$(wc -c < "$meta_temp_dir/sa_rate_index/test_ref.bwt.2bit.64")
log "Index size: default (-u 3) = $default_size bytes, -u 1 = $dense_size bytes"

if [[ "$dense_size" -le "$default_size" ]]; then
  log_error "Expected -u 1 index to be larger than default (-u 3) index"
  exit 1
fi

log "✅ TEST 5 completed successfully"

print_test_summary "All tests completed successfully"

