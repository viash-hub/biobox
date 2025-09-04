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

# --- Test Case 1: Basic intersection ---
log "Starting TEST 1: Basic intersection"

# Create test BED files
log "Creating test BED data..."
cat > "$test_dir/featuresA.bed" << 'EOF'
chr1	100	200	feature1
chr1	300	400	feature2
chr2	500	600	feature3
EOF

cat > "$test_dir/featuresB.bed" << 'EOF'
chr1	150	250	overlapping1
chr1	350	450	overlapping2
chr2	550	650	overlapping3
EOF

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input_a "$test_dir/featuresA.bed" \
  --input_b "$test_dir/featuresB.bed" \
  --output "$meta_temp_dir/output1.bed"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output1.bed" "output intersection file"
check_file_not_empty "$meta_temp_dir/output1.bed" "output intersection file"
check_file_contains "$meta_temp_dir/output1.bed" "chr1"
log "✅ TEST 1 completed successfully"

# --- Test Case 2: Intersection with -wa option ---
log "Starting TEST 2: Intersection with -wa (write A) option"

log "Executing $meta_name with -wa option..."
"$meta_executable" \
  --input_a "$test_dir/featuresA.bed" \
  --input_b "$test_dir/featuresB.bed" \
  --write_a \
  --output "$meta_temp_dir/output2.bed"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output2.bed" "output file with -wa"
check_file_not_empty "$meta_temp_dir/output2.bed" "output file with -wa"
check_file_contains "$meta_temp_dir/output2.bed" "feature"
log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intersection with -wb option ---
log "Starting TEST 3: Intersection with -wb (write B) option"

log "Executing $meta_name with -wb option..."
"$meta_executable" \
  --input_a "$test_dir/featuresA.bed" \
  --input_b "$test_dir/featuresB.bed" \
  --write_b \
  --output "$meta_temp_dir/output3.bed"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output3.bed" "output file with -wb"
check_file_not_empty "$meta_temp_dir/output3.bed" "output file with -wb"
check_file_contains "$meta_temp_dir/output3.bed" "overlapping"
log "✅ TEST 3 completed successfully"
