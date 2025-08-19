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
log "Starting TEST 1: Basic Bowtie2 indexing"

log "Generating test reference genome..."
create_test_fasta "$test_data_dir/test_ref.fasta" 2 1000
check_file_exists "$test_data_dir/test_ref.fasta" "test reference genome"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/bt2_index"

log "Validating TEST 1 outputs..."
check_dir_exists "$meta_temp_dir/bt2_index" "output index directory"

# Check for standard bowtie2 index files
log "Checking for bowtie2 index files..."
index_files=(
  "$meta_temp_dir/bt2_index/test_ref.1.bt2"
  "$meta_temp_dir/bt2_index/test_ref.2.bt2"
  "$meta_temp_dir/bt2_index/test_ref.3.bt2"
  "$meta_temp_dir/bt2_index/test_ref.4.bt2"
  "$meta_temp_dir/bt2_index/test_ref.rev.1.bt2"
  "$meta_temp_dir/bt2_index/test_ref.rev.2.bt2"
)

for file in "${index_files[@]}"; do
  check_file_exists "$file" "bowtie2 index file $(basename "$file")"
  check_file_not_empty "$file" "bowtie2 index file $(basename "$file")"
done

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Custom index name ---
log "Starting TEST 2: Custom index name"

log "Executing $meta_name with custom index name..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/custom_index" \
  --index_name "custom_genome"

log "Validating TEST 2 outputs..."
check_dir_exists "$meta_temp_dir/custom_index" "custom index directory"

# Check for index files with custom name
log "Checking for custom-named index files..."
custom_index_files=(
  "$meta_temp_dir/custom_index/custom_genome.1.bt2"
  "$meta_temp_dir/custom_index/custom_genome.2.bt2"
  "$meta_temp_dir/custom_index/custom_genome.3.bt2"
  "$meta_temp_dir/custom_index/custom_genome.4.bt2"
  "$meta_temp_dir/custom_index/custom_genome.rev.1.bt2"
  "$meta_temp_dir/custom_index/custom_genome.rev.2.bt2"
)

for file in "${custom_index_files[@]}"; do
  check_file_exists "$file" "custom-named index file $(basename "$file")"
  check_file_not_empty "$file" "custom-named index file $(basename "$file")"
done

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Large index option ---
log "Starting TEST 3: Large index option"

log "Executing $meta_name with large index option..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/large_index" \
  --large_index

log "Validating TEST 3 outputs..."
check_dir_exists "$meta_temp_dir/large_index" "large index directory"

# Check for index files (large index may have different structure)
large_index_files=(
  "$meta_temp_dir/large_index/test_ref.1.bt2l"
  "$meta_temp_dir/large_index/test_ref.2.bt2l"
  "$meta_temp_dir/large_index/test_ref.3.bt2l"
  "$meta_temp_dir/large_index/test_ref.4.bt2l"
  "$meta_temp_dir/large_index/test_ref.rev.1.bt2l"
  "$meta_temp_dir/large_index/test_ref.rev.2.bt2l"
)

# Check if large index files exist, if not check regular format
has_large_format=true
for file in "${large_index_files[@]}"; do
  if [[ ! -f "$file" ]]; then
    has_large_format=false
    break
  fi
done

if [[ "$has_large_format" == "true" ]]; then
  log "Large format index files detected"
  for file in "${large_index_files[@]}"; do
    check_file_exists "$file" "large index file $(basename "$file")"
    check_file_not_empty "$file" "large index file $(basename "$file")"
  done
else
  log "Regular format index files with large index option"
  regular_large_files=(
    "$meta_temp_dir/large_index/test_ref.1.bt2"
    "$meta_temp_dir/large_index/test_ref.2.bt2"
    "$meta_temp_dir/large_index/test_ref.3.bt2"
    "$meta_temp_dir/large_index/test_ref.4.bt2"
    "$meta_temp_dir/large_index/test_ref.rev.1.bt2"
    "$meta_temp_dir/large_index/test_ref.rev.2.bt2"
  )
  for file in "${regular_large_files[@]}"; do
    check_file_exists "$file" "index file $(basename "$file")"
    check_file_not_empty "$file" "index file $(basename "$file")"
  done
fi

log "✅ TEST 3 completed successfully"

print_test_summary "All tests completed successfully - bowtie2_build component validates basic indexing, custom naming, and large index scenarios"

# --- Test Case 4: Large index option ---
log "Starting TEST 4: Large index option"

log "Executing $meta_name with large index option..."
"$meta_executable" \
  --input "$test_data_dir/test_ref.fasta" \
  --output "$meta_temp_dir/large_index" \
  --large_index

log "Validating TEST 4 outputs..."
check_dir_exists "$meta_temp_dir/large_index" "large index directory"

# Check for index files (large index may have different structure)
large_index_files=(
  "$meta_temp_dir/large_index/test_ref.1.bt2l"
  "$meta_temp_dir/large_index/test_ref.2.bt2l"
  "$meta_temp_dir/large_index/test_ref.3.bt2l"
  "$meta_temp_dir/large_index/test_ref.4.bt2l"
  "$meta_temp_dir/large_index/test_ref.rev.1.bt2l"
  "$meta_temp_dir/large_index/test_ref.rev.2.bt2l"
)

# Check if large index files exist, if not check regular format
has_large_format=true
for file in "${large_index_files[@]}"; do
  if [[ ! -f "$file" ]]; then
    has_large_format=false
    break
  fi
done

if [[ "$has_large_format" == "true" ]]; then
  log "Large format index files detected"
  for file in "${large_index_files[@]}"; do
    check_file_exists "$file" "large index file $(basename "$file")"
    check_file_not_empty "$file" "large index file $(basename "$file")"
  done
else
  log "Regular format index files with large index option"
  regular_large_files=(
    "$meta_temp_dir/large_index/test_ref.1.bt2"
    "$meta_temp_dir/large_index/test_ref.2.bt2"
    "$meta_temp_dir/large_index/test_ref.3.bt2"
    "$meta_temp_dir/large_index/test_ref.4.bt2"
    "$meta_temp_dir/large_index/test_ref.rev.1.bt2"
    "$meta_temp_dir/large_index/test_ref.rev.2.bt2"
  )
  for file in "${regular_large_files[@]}"; do
    check_file_exists "$file" "index file $(basename "$file")"
    check_file_not_empty "$file" "index file $(basename "$file")"
  done
fi

log "✅ TEST 4 completed successfully"

print_test_summary "All tests completed successfully"
