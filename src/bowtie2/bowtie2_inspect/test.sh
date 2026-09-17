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

# Prepare test data
log "Generating test reference genome..."
create_test_fasta "$test_data_dir/test_ref.fasta" 2 1000
check_file_exists "$test_data_dir/test_ref.fasta" "test reference genome"

# Build index using bowtie2-build
log "Building Bowtie2 index for inspection tests..."
mkdir -p "$test_data_dir/index"
bowtie2-build "$test_data_dir/test_ref.fasta" "$test_data_dir/index/test_ref" >/dev/null 2>&1

# Verify index was created
check_file_exists "$test_data_dir/index/test_ref.1.bt2" "bowtie2 index file"

# --- Test Case 1: Default FASTA output ---
log "Starting TEST 1: Default FASTA output"

log "Executing $meta_name with default FASTA output..."
"$meta_executable" \
  --index "$test_data_dir/index" \
  --output "$meta_temp_dir/sequences.fasta"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/sequences.fasta" "FASTA output"
check_file_not_empty "$meta_temp_dir/sequences.fasta" "FASTA output"

# Check FASTA format
if grep -q "^>" "$meta_temp_dir/sequences.fasta"; then
  log "✓ Output contains FASTA headers"
else
  log_error "Output does not contain proper FASTA headers"
  exit 1
fi

# Check for sequence content
if grep -q "^[ATCGN]" "$meta_temp_dir/sequences.fasta"; then
  log "✓ Output contains nucleotide sequences"
else
  log_error "Output does not contain nucleotide sequences"
  exit 1
fi

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Names only output ---
log "Starting TEST 2: Names only output"

log "Executing $meta_name with names only..."
"$meta_executable" \
  --index "$test_data_dir/index" \
  --names \
  --output "$meta_temp_dir/names.txt"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/names.txt" "names output"
check_file_not_empty "$meta_temp_dir/names.txt" "names output"

# Check that output contains sequence names from our test FASTA
if grep -q "seq" "$meta_temp_dir/names.txt"; then
  log "✓ Output contains expected sequence names"
else
  log_error "Output does not contain expected sequence names"
  exit 1
fi

# Ensure it doesn't contain sequence data (should be names only)
if ! grep -q "^[ATCGN]" "$meta_temp_dir/names.txt"; then
  log "✓ Output correctly contains only names, no sequences"
else
  log_error "Output incorrectly contains sequence data"
  exit 1
fi

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Summary output ---
log "Starting TEST 3: Summary output"

log "Executing $meta_name with summary..."
"$meta_executable" \
  --index "$test_data_dir/index" \
  --summary \
  --output "$meta_temp_dir/summary.txt"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/summary.txt" "summary output"
check_file_not_empty "$meta_temp_dir/summary.txt" "summary output"

# Check for summary-specific content
if grep -q -i "sequence\|length\|total" "$meta_temp_dir/summary.txt"; then
  log "✓ Output contains summary information"
else
  log_error "Output does not contain expected summary information"
  exit 1
fi

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Standard output (no output file) ---
log "Starting TEST 4: Standard output"

log "Executing $meta_name with stdout output..."
stdout_output=$("$meta_executable" --index "$test_data_dir/index" --names 2>/dev/null)

log "Validating TEST 4 outputs..."
if [[ -n "$stdout_output" ]]; then
  log "✓ Standard output contains data"
else
  log_error "Standard output is empty"
  exit 1
fi

# Check that stdout contains expected content
if echo "$stdout_output" | grep -q "seq"; then
  log "✓ Standard output contains expected sequence names"
else
  log_error "Standard output does not contain expected content"
  exit 1
fi

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Across parameter ---
log "Starting TEST 5: Across parameter"

log "Executing $meta_name with across parameter..."
"$meta_executable" \
  --index "$test_data_dir/index" \
  --across 60 \
  --output "$meta_temp_dir/across.fasta"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/across.fasta" "across output"
check_file_not_empty "$meta_temp_dir/across.fasta" "across output"

# Check FASTA format
if grep -q "^>" "$meta_temp_dir/across.fasta"; then
  log "✓ Across output contains FASTA headers"
else
  log_error "Across output does not contain proper FASTA headers"
  exit 1
fi

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Large index directory ---
log "Starting TEST 6: Large index directory"

log "Building a 'large' bowtie2 index (.bt2l files)..."
mkdir -p "$test_data_dir/large_index"
bowtie2-build --large-index "$test_data_dir/test_ref.fasta" "$test_data_dir/large_index/test_ref" >/dev/null 2>&1
check_file_exists "$test_data_dir/large_index/test_ref.1.bt2l" "large bowtie2 index file"

log "Executing $meta_name with a large index directory..."
"$meta_executable" \
  --index "$test_data_dir/large_index" \
  --names \
  --output "$meta_temp_dir/large_index_names.txt"

log "Validating TEST 6 outputs..."
check_file_exists "$meta_temp_dir/large_index_names.txt" "large index names output"
check_file_not_empty "$meta_temp_dir/large_index_names.txt" "large index names output"
check_file_contains "$meta_temp_dir/large_index_names.txt" "seq" "large index names output"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: Directory without index files ---
log "Starting TEST 7: Directory without bowtie2 index files"

mkdir -p "$test_data_dir/empty_index"

log "Executing $meta_name with an index directory that holds no index files..."
if "$meta_executable" \
  --index "$test_data_dir/empty_index" \
  --names \
  --output "$meta_temp_dir/empty_index.txt" 2> "$meta_temp_dir/empty_index.log"; then
  log_error "Expected failure when the index directory contains no index files"
  exit 1
fi

log "Validating TEST 7 outputs..."
check_file_contains "$meta_temp_dir/empty_index.log" "no bowtie2 index files" "error message about missing index files"

log "✅ TEST 7 completed successfully"

# --- Test Case 8: Directory with multiple indices ---
log "Starting TEST 8: Directory containing multiple bowtie2 indices"

mkdir -p "$test_data_dir/multi_index"
for file in "$test_data_dir/index/test_ref".*.bt2; do
  cp "$file" "$test_data_dir/multi_index/"
  cp "$file" "$test_data_dir/multi_index/other$(basename "$file" | sed 's/^test_ref//')"
done

log "Executing $meta_name with an index directory that holds two indices..."
if "$meta_executable" \
  --index "$test_data_dir/multi_index" \
  --names \
  --output "$meta_temp_dir/multi_index.txt" 2> "$meta_temp_dir/multi_index.log"; then
  log_error "Expected failure when the index directory contains multiple indices"
  exit 1
fi

log "Validating TEST 8 outputs..."
check_file_contains "$meta_temp_dir/multi_index.log" "multiple bowtie2 indices" "error message about multiple indices"
check_file_contains "$meta_temp_dir/multi_index.log" "index_prefix" "error message pointing at --index_prefix"

log "✅ TEST 8 completed successfully"

# --- Test Case 9: Index directory of symlinks ---
log "Starting TEST 9: Index directory containing symlinks"

log "Creating an index directory of symlinks, as a workflow engine would stage it..."
mkdir -p "$test_data_dir/symlink_index"
ln -s "$test_data_dir/index/test_ref".*.bt2 "$test_data_dir/symlink_index/"

log "Executing $meta_name with a symlinked index directory..."
"$meta_executable" \
  --index "$test_data_dir/symlink_index" \
  --names \
  --output "$meta_temp_dir/symlink_index_names.txt"

log "Validating TEST 9 outputs..."
check_file_exists "$meta_temp_dir/symlink_index_names.txt" "symlinked index names output"
check_file_contains "$meta_temp_dir/symlink_index_names.txt" "seq" "symlinked index names output"

log "✅ TEST 9 completed successfully"

# --- Test Case 10: Explicit --large_index with a large index ---
log "Starting TEST 10: Explicit --large_index with a large index"

log "Executing $meta_name with --large_index on a large index directory..."
"$meta_executable" \
  --index "$test_data_dir/large_index" \
  --large_index \
  --names \
  --output "$meta_temp_dir/forced_large_names.txt"

log "Validating TEST 10 outputs..."
check_file_exists "$meta_temp_dir/forced_large_names.txt" "forced large index names output"
check_file_contains "$meta_temp_dir/forced_large_names.txt" "seq" "forced large index names output"

log "✅ TEST 10 completed successfully"

# --- Test Case 11: Explicit --large_index without a large index ---
log "Starting TEST 11: Explicit --large_index without a large index"

log "Executing $meta_name with --large_index on a small index directory..."
if "$meta_executable" \
  --index "$test_data_dir/index" \
  --large_index \
  --names \
  --output "$meta_temp_dir/forced_large_missing.txt" 2> "$meta_temp_dir/forced_large_missing.log"; then
  log_error "Expected failure when --large_index is set but no .bt2l files are present"
  exit 1
fi

log "Validating TEST 11 outputs..."
check_file_contains "$meta_temp_dir/forced_large_missing.log" "no large bowtie2 index files" \
  "error message about missing large index files"

log "✅ TEST 11 completed successfully"

# --- Test Case 12: Selecting one index out of a directory holding several ---
log "Starting TEST 12: Selecting an index with --index_prefix"

# The directory built for TEST 8 holds the 'test_ref' and 'other' indices, so the
# prefix is what decides which of the two is inspected.
log "Executing $meta_name with --index_prefix on an ambiguous index directory..."
"$meta_executable" \
  --index "$test_data_dir/multi_index" \
  --index_prefix "other" \
  --names \
  --output "$meta_temp_dir/index_prefix_names.txt"

log "Validating TEST 12 outputs..."
check_file_exists "$meta_temp_dir/index_prefix_names.txt" "selected-index names output"
check_file_not_empty "$meta_temp_dir/index_prefix_names.txt" "selected-index names output"
check_file_contains "$meta_temp_dir/index_prefix_names.txt" "seq" "selected-index names output"

log "✅ TEST 12 completed successfully"

# --- Test Case 13: An index prefix that is not present ---
log "Starting TEST 13: Index prefix that is not present"

log "Executing $meta_name with an index prefix that is not present..."
if "$meta_executable" \
  --index "$test_data_dir/index" \
  --index_prefix "not_an_index" \
  --names \
  --output "$meta_temp_dir/unknown_prefix.txt" 2> "$meta_temp_dir/unknown_prefix.log"; then
  log_error "Expected failure when the index prefix is not present in the index directory"
  exit 1
fi

log "Validating TEST 13 outputs..."
check_file_contains "$meta_temp_dir/unknown_prefix.log" "not_an_index" "error message naming the index prefix"

log "✅ TEST 13 completed successfully"

# --- Test Case 14: --index_prefix combined with --large_index ---
log "Starting TEST 14: --index_prefix combined with --large_index"

log "Building a second large index next to the first one..."
bowtie2-build --large-index "$test_data_dir/test_ref.fasta" "$test_data_dir/large_index/other" >/dev/null 2>&1
check_file_exists "$test_data_dir/large_index/other.1.bt2l" "second large bowtie2 index file"

log "Executing $meta_name with --large_index and --index_prefix..."
"$meta_executable" \
  --index "$test_data_dir/large_index" \
  --large_index \
  --index_prefix "other" \
  --names \
  --output "$meta_temp_dir/large_index_prefix_names.txt"

log "Validating TEST 14 outputs..."
check_file_exists "$meta_temp_dir/large_index_prefix_names.txt" "selected large index names output"
check_file_contains "$meta_temp_dir/large_index_prefix_names.txt" "seq" "selected large index names output"

log "✅ TEST 14 completed successfully"

print_test_summary "All tests completed successfully"
