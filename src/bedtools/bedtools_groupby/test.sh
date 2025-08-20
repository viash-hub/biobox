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

# Create test directory
test_dir="$meta_temp_dir/test_data"
mkdir -p "$test_dir"

# Create test BED file with data for grouping
log "Creating test BED data..."
cat > "$test_dir/test.bed" << 'EOF'
chr1	100	200	feature1	10	+
chr1	300	400	feature2	20	+
chr1	500	600	feature3	30	+
chr2	100	200	feature4	15	-
chr2	300	400	feature5	25	-
chr3	100	200	feature6	35	+
EOF

# --- Test Case 1: Basic grouping by column 1 (chromosome) with sum operation ---
log "Starting TEST 1: Basic grouping by chromosome with sum"

log "Executing $meta_name with basic grouping..."
"$meta_executable" \
  --input "$test_dir/test.bed" \
  --groupby 1 \
  --column 5 \
  --operation sum \
  --output "$meta_temp_dir/output1.txt"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output1.txt" "grouped output file"
check_file_not_empty "$meta_temp_dir/output1.txt" "grouped output file"
check_file_contains "$meta_temp_dir/output1.txt" "chr1"
check_file_contains "$meta_temp_dir/output1.txt" "chr2"
check_file_contains "$meta_temp_dir/output1.txt" "chr3"
log "✅ TEST 1 completed successfully"

# --- Test Case 2: Group by multiple columns with mean operation ---
log "Starting TEST 2: Group by chromosome and strand with mean"

log "Executing $meta_name with multiple column grouping..."
"$meta_executable" \
  --input "$test_dir/test.bed" \
  --groupby 1,6 \
  --column 5 \
  --operation mean \
  --output "$meta_temp_dir/output2.txt"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output2.txt" "multi-column grouped output"
check_file_not_empty "$meta_temp_dir/output2.txt" "multi-column grouped output"
check_file_contains "$meta_temp_dir/output2.txt" "chr1"
check_file_contains "$meta_temp_dir/output2.txt" "+"
check_file_contains "$meta_temp_dir/output2.txt" "-"
log "✅ TEST 2 completed successfully"

# --- Test Case 3: Count operation ---
log "Starting TEST 3: Group by chromosome with count operation"

log "Executing $meta_name with count operation..."
"$meta_executable" \
  --input "$test_dir/test.bed" \
  --groupby 1 \
  --column 5 \
  --operation count \
  --output "$meta_temp_dir/output3.txt"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output3.txt" "count output file"
check_file_not_empty "$meta_temp_dir/output3.txt" "count output file"
# chr1 should have 3 features, chr2 should have 2, chr3 should have 1
check_file_contains "$meta_temp_dir/output3.txt" "3"
check_file_contains "$meta_temp_dir/output3.txt" "2"
check_file_contains "$meta_temp_dir/output3.txt" "1"
log "✅ TEST 3 completed successfully"

# --- Test Case 4: Min/Max operations ---
log "Starting TEST 4: Group by chromosome with min operation"

log "Executing $meta_name with min operation..."
"$meta_executable" \
  --input "$test_dir/test.bed" \
  --groupby 1 \
  --column 5 \
  --operation min \
  --output "$meta_temp_dir/output4.txt"

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output4.txt" "min output file"
check_file_not_empty "$meta_temp_dir/output4.txt" "min output file"
log "✅ TEST 4 completed successfully"

# --- Test Case 5: Full output with additional options ---
log "Starting TEST 5: Group with full output and header"

log "Executing $meta_name with full output options..."
"$meta_executable" \
  --input "$test_dir/test.bed" \
  --groupby 1 \
  --column 5 \
  --operation sum \
  --full \
  --output "$meta_temp_dir/output5.txt"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/output5.txt" "full output file"
check_file_not_empty "$meta_temp_dir/output5.txt" "full output file"
# Full output should include more columns from original data
log "✅ TEST 5 completed successfully"

log "🎉 All tests completed successfully for $meta_name!"
