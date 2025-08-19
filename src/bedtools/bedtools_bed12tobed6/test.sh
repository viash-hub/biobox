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

# Create a test BED12 file
log "Creating test BED12 data..."
cat > "$test_dir/test.bed12" << 'EOF'
chr1	100	600	gene1	1000	+	100	600	255,0,0	3	100,150,200	0,200,300
chr2	200	800	gene2	800	-	200	800	0,255,0	2	200,250	0,350
chr3	300	500	gene3	500	.	300	500	0,0,255	1	200	0
EOF

# --- Test Case 1: Basic BED12 to BED6 conversion ---
log "Starting TEST 1: Basic BED12 to BED6 conversion"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
    --input "$test_dir/test.bed12" \
    --output "$meta_temp_dir/output1.bed6"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output1.bed6" "output BED6 file"
check_file_not_empty "$meta_temp_dir/output1.bed6" "output BED6 file"

# Check that BED6 file has correct number of columns (6 columns)
awk 'NF != 6 { exit 1 }' "$meta_temp_dir/output1.bed6" || { 
    log_error "Output is not in BED6 format (expected 6 columns per line)"
    exit 1 
}

# Check that we have more BED6 entries than BED12 entries (due to block splitting)
bed12_lines=$(wc -l < "$test_dir/test.bed12")
bed6_lines=$(wc -l < "$meta_temp_dir/output1.bed6")
log "Input BED12: $bed12_lines lines, Output BED6: $bed6_lines lines"

[ "$bed6_lines" -gt "$bed12_lines" ] || { 
    log_error "Expected more BED6 lines than BED12 lines due to block splitting"
    exit 1 
}

# Check that gene names are preserved
check_file_contains "$meta_temp_dir/output1.bed6" "gene1" "gene names from BED12"
check_file_contains "$meta_temp_dir/output1.bed6" "gene2" "gene names from BED12"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: BED12 to BED6 with --n_score option ---
log "Starting TEST 2: BED12 to BED6 with block numbering"

log "Executing $meta_name with --n_score flag..."
"$meta_executable" \
    --input "$test_dir/test.bed12" \
    --output "$meta_temp_dir/output2.bed6" \
    --n_score

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output2.bed6" "output BED6 file with block numbers"
check_file_not_empty "$meta_temp_dir/output2.bed6" "output BED6 file with block numbers"

# Check that BED6 file has correct number of columns
awk 'NF != 6 { exit 1 }' "$meta_temp_dir/output2.bed6" || { 
    log_error "Output is not in BED6 format (expected 6 columns per line)"
    exit 1 
}

# Check that scores are block numbers (should contain "1", "2", "3" for gene1 with 3 blocks)
check_file_contains "$meta_temp_dir/output2.bed6" $'\t1\t' "block number 1 in score column"
check_file_contains "$meta_temp_dir/output2.bed6" $'\t2\t' "block number 2 in score column"
check_file_contains "$meta_temp_dir/output2.bed6" $'\t3\t' "block number 3 in score column"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Test with single-block BED12 ---
log "Starting TEST 3: Single-block BED12 conversion"

# Create a simple single-block BED12 (should produce single BED6)
cat > "$test_dir/single_block.bed12" << 'EOF'
chrX	1000	2000	single_gene	900	+	1000	2000	128,128,128	1	1000	0
EOF

log "Executing $meta_name with single-block BED12..."
"$meta_executable" \
    --input "$test_dir/single_block.bed12" \
    --output "$meta_temp_dir/output3.bed6"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output3.bed6" "single-block BED6 output"
check_file_not_empty "$meta_temp_dir/output3.bed6" "single-block BED6 output"

# Should have exactly one line (single block)
single_lines=$(wc -l < "$meta_temp_dir/output3.bed6")
[ "$single_lines" -eq 1 ] || { 
    log_error "Expected exactly 1 line for single-block BED12, got $single_lines"
    exit 1 
}

# Check that it contains the expected gene name
check_file_contains "$meta_temp_dir/output3.bed6" "single_gene" "single gene name"

log "✅ TEST 3 completed successfully"

print_test_summary "All tests completed successfully"
