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

# Create test FASTA file
log "Creating test FASTA data..."
cat > "$test_dir/test.fa" << 'EOF'
>chr1
AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG
>chr2
TTTTTTTTGGGGGGGGGGGGGGCGGATCGGGGGGGGGGGGGGAAA
EOF

# Create test BED file  
cat > "$test_dir/test.bed" << 'EOF'
chr1	5	10	seq1
chr2	15	20	seq2
EOF

# --- Test Case 1: Basic FASTA sequence extraction ---
log "Starting TEST 1: Basic FASTA sequence extraction"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
    --input_bed "$test_dir/test.bed" \
    --input_fasta "$test_dir/test.fa" \
    --output "$meta_temp_dir/output1.fasta"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output1.fasta" "output FASTA file"
check_file_not_empty "$meta_temp_dir/output1.fasta" "output FASTA file"
check_file_contains "$meta_temp_dir/output1.fasta" ">chr1:5-10"
check_file_contains "$meta_temp_dir/output1.fasta" "AAACC"
log "✅ TEST 1 completed successfully"

# --- Test Case 2: FASTA extraction with --name option ---
log "Starting TEST 2: FASTA extraction with --name option"

log "Executing $meta_name with --name option..."
"$meta_executable" \
    --input_bed "$test_dir/test.bed" \
    --input_fasta "$test_dir/test.fa" \
    --name \
    --output "$meta_temp_dir/output2.fasta"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output2.fasta" "output FASTA file with names"
check_file_not_empty "$meta_temp_dir/output2.fasta" "output FASTA file with names"
check_file_contains "$meta_temp_dir/output2.fasta" ">seq1::chr1:5-10"
check_file_contains "$meta_temp_dir/output2.fasta" ">seq2::chr2:15-20"
log "✅ TEST 2 completed successfully"

# --- Test Case 3: FASTA extraction with --name_only option ---
log "Starting TEST 3: FASTA extraction with --name_only option"

log "Executing $meta_name with --name_only option..."
"$meta_executable" \
    --input_bed "$test_dir/test.bed" \
    --input_fasta "$test_dir/test.fa" \
    --name_only \
    --output "$meta_temp_dir/output3.fasta"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output3.fasta" "output FASTA file with name only"
check_file_not_empty "$meta_temp_dir/output3.fasta" "output FASTA file with name only"
check_file_contains "$meta_temp_dir/output3.fasta" ">seq1"
check_file_contains "$meta_temp_dir/output3.fasta" ">seq2"
log "✅ TEST 3 completed successfully"

# --- Test Case 4: Tab-delimited output ---
log "Starting TEST 4: Tab-delimited output with --tab option"

log "Executing $meta_name with --tab option..."
"$meta_executable" \
    --input_bed "$test_dir/test.bed" \
    --input_fasta "$test_dir/test.fa" \
    --name_only \
    --tab \
    --output "$meta_temp_dir/output4.txt"

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output4.txt" "tab-delimited output file"
check_file_not_empty "$meta_temp_dir/output4.txt" "tab-delimited output file"
check_file_contains "$meta_temp_dir/output4.txt" "seq1"
check_file_contains "$meta_temp_dir/output4.txt" "AAACC"
log "✅ TEST 4 completed successfully"

# --- Test Case 5: BED output format ---
log "Starting TEST 5: BED output format with --bed_out option"

log "Executing $meta_name with --bed_out option..."
"$meta_executable" \
    --input_bed "$test_dir/test.bed" \
    --input_fasta "$test_dir/test.fa" \
    --bed_out \
    --output "$meta_temp_dir/output5.bed"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/output5.bed" "BED output file"
check_file_not_empty "$meta_temp_dir/output5.bed" "BED output file"
# BED format output contains sequences with coordinates
log "✅ TEST 5 completed successfully"

log "🎉 All tests completed successfully for $meta_name!"
