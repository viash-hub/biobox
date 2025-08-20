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

# Create a test SAM file with proper format (based on original test data)
log "Creating test SAM data..."
cat > "$test_dir/test.sam" << 'EOF'
@SQ	SN:chr2:172936693-172938111	LN:1418
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188
my_read/1	99	chr2:172936693-172938111	129	60	100M	=	429	400	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	SM:i:85
my_read/2	147	chr2:172936693-172938111	429	60	100M	=	129	-400	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	SM:i:85
EOF

# Convert SAM to BAM using samtools (if available in container) or use the SAM directly
log "Converting SAM to BAM..."
if command -v samtools >/dev/null 2>&1; then
  samtools view -bS "$test_dir/test.sam" > "$test_dir/test.bam"
  input_file="$test_dir/test.bam"
else
  # bedtools can handle SAM files directly
  input_file="$test_dir/test.sam"
  log "Using SAM file directly (samtools not available)"
fi

# --- Test Case 1: Basic BAM to BED conversion ---
log "Starting TEST 1: Basic BAM to BED conversion"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$input_file" \
  --output "$meta_temp_dir/output1.bed"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output1.bed" "output BED file"
check_file_not_empty "$meta_temp_dir/output1.bed" "output BED file"

# Check that BED file has correct number of columns (6 for BED6)
line_count=$(wc -l < "$meta_temp_dir/output1.bed")
log "Output contains $line_count lines"
[ "$line_count" -gt 0 ] || { log_error "Output file is empty"; exit 1; }

# Check that each line has 6 columns (BED6 format)
awk 'NF != 6 { exit 1 }' "$meta_temp_dir/output1.bed" || { 
  log_error "Output is not in BED6 format (expected 6 columns per line)"
  exit 1 
}

log "✅ TEST 1 completed successfully"

# --- Test Case 2: BEDPE format ---
log "Starting TEST 2: BEDPE format conversion"

log "Executing $meta_name with --bedpe flag..."
"$meta_executable" \
  --input "$input_file" \
  --output "$meta_temp_dir/output2.bedpe" \
  --bedpe

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output2.bedpe" "output BEDPE file"
check_file_not_empty "$meta_temp_dir/output2.bedpe" "output BEDPE file"

# Check that BEDPE file has correct number of columns (10 for BEDPE)
awk 'NF != 10 { exit 1 }' "$meta_temp_dir/output2.bedpe" || { 
  log_error "Output is not in BEDPE format (expected 10 columns per line)"
  exit 1 
}

log "✅ TEST 2 completed successfully"

# --- Test Case 3: BED12 format ---
log "Starting TEST 3: BED12 format conversion"

log "Executing $meta_name with --bed12 flag..."
"$meta_executable" \
  --input "$input_file" \
  --output "$meta_temp_dir/output3.bed12" \
  --bed12

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output3.bed12" "output BED12 file"
check_file_not_empty "$meta_temp_dir/output3.bed12" "output BED12 file"

# Check that BED12 file has correct number of columns (12 for BED12)
awk 'NF != 12 { exit 1 }' "$meta_temp_dir/output3.bed12" || { 
  log_error "Output is not in BED12 format (expected 12 columns per line)"
  exit 1 
}

log "✅ TEST 3 completed successfully"

# --- Test Case 4: CIGAR addition ---
log "Starting TEST 4: CIGAR string addition"

log "Executing $meta_name with --cigar flag..."
"$meta_executable" \
  --input "$input_file" \
  --output "$meta_temp_dir/output4.bed" \
  --cigar

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output4.bed" "output BED file with CIGAR"
check_file_not_empty "$meta_temp_dir/output4.bed" "output BED file with CIGAR"

# Check that BED file has correct number of columns (7 for BED6 + CIGAR)
awk 'NF != 7 { exit 1 }' "$meta_temp_dir/output4.bed" || { 
  log_error "Output is not in BED6+CIGAR format (expected 7 columns per line)"
  exit 1 
}

# Check that the 7th column contains CIGAR strings
check_file_contains "$meta_temp_dir/output4.bed" "100M" "BED file with CIGAR strings"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests completed successfully"
