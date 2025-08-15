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
  --index "$test_data_dir/index/test_ref" \
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
  --index "$test_data_dir/index/test_ref" \
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
  --index "$test_data_dir/index/test_ref" \
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
stdout_output=$("$meta_executable" --index "$test_data_dir/index/test_ref" --names 2>/dev/null)

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
  --index "$test_data_dir/index/test_ref" \
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

print_test_summary "All tests completed successfully"
