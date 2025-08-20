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

# Create a test SAM file with proper FASTQ data
log "Creating test SAM data..."
cat > "$test_dir/test.sam" << 'EOF'
@SQ	SN:chr1	LN:1000
@PG	ID:bwa	PN:bwa	VN:0.7.17
read1	0	chr1	100	60	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	0	chr1	200	60	50M	*	0	0	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
EOF

# --- Test Case 1: Basic BAM to FASTQ conversion (single-end) ---
log "Starting TEST 1: Basic BAM to FASTQ conversion"

log "Executing $meta_name with single-end BAM..."
"$meta_executable" \
  --input "$test_dir/test.sam" \
  --fastq "$meta_temp_dir/output1.fastq"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output1.fastq" "output FASTQ file"
check_file_not_empty "$meta_temp_dir/output1.fastq" "output FASTQ file"

# Check FASTQ format (should have 4 lines per read: header, sequence, +, quality)
total_lines=$(wc -l < "$meta_temp_dir/output1.fastq")
log "Output FASTQ contains $total_lines lines"
[ $((total_lines % 4)) -eq 0 ] || { log_error "FASTQ format error: line count not divisible by 4"; exit 1; }

# Check that FASTQ contains expected patterns
check_file_contains "$meta_temp_dir/output1.fastq" "@read1" "FASTQ headers"
check_file_contains "$meta_temp_dir/output1.fastq" "AAAAAAAA" "sequence content"
check_file_contains "$meta_temp_dir/output1.fastq" "IIIIIIII" "quality scores"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Test --tags option ---
log "Starting TEST 2: BAM to FASTQ with --tags option"

# For the tags test, we'll just verify the command runs without error
# since creating BAM with R2/Q2 tags would be complex
log "Executing $meta_name with --tags flag..."
"$meta_executable" \
  --input "$test_dir/test.sam" \
  --fastq "$meta_temp_dir/output2.fastq" \
  --tags

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output2.fastq" "output FASTQ file with tags"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Test with secondary output (without actual paired data) ---
log "Starting TEST 3: Test secondary output parameter"

# Test that the fastq2 parameter is accepted (even if no paired reads are present)
log "Executing $meta_name with --fastq2 parameter..."
"$meta_executable" \
  --input "$test_dir/test.sam" \
  --fastq "$meta_temp_dir/output3_R1.fastq" \
  --fastq2 "$meta_temp_dir/output3_R2.fastq"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output3_R1.fastq" "primary FASTQ file"
check_file_not_empty "$meta_temp_dir/output3_R1.fastq" "primary FASTQ file"

# The R2 file may be empty since we don't have paired reads, but should exist
check_file_exists "$meta_temp_dir/output3_R2.fastq" "secondary FASTQ file"

log "✅ TEST 3 completed successfully"

print_test_summary "All tests completed successfully"


