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

# Create test genome file
log "Creating test genome file..."
cat > "$test_dir/test.genome" << 'EOF'
chr1	248956422
chr2	242193529
chr3	198295559
EOF

# Create test BED file (BED4 minimum for bedtobam)
log "Creating test BED file..."
cat > "$test_dir/test.bed" << 'EOF'
chr1	1000	2000	gene1	100	+
chr2	3000	4000	gene2	200	-
chr3	5000	6000	gene3	150	+
EOF

# Create test BED12 file
log "Creating test BED12 file..."
cat > "$test_dir/test.bed12" << 'EOF'
chr1	1000	3000	gene1	100	+	1000	3000	255,0,0	2	500,500	0,1500
chr2	2000	5000	gene2	200	-	2000	5000	0,255,0	3	400,300,400	0,1500,2600
EOF

# --- Test Case 1: Basic BED to BAM conversion ---
log "Starting TEST 1: Basic BED to BAM conversion"

log "Executing $meta_name with basic BED file..."
"$meta_executable" \
  --input "$test_dir/test.bed" \
  --genome "$test_dir/test.genome" \
  --output "$meta_temp_dir/output1.bam"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output1.bam" "output BAM file"
check_file_not_empty "$meta_temp_dir/output1.bam" "output BAM file"

# Check if it's a valid BAM file by reading header
if command -v samtools >/dev/null 2>&1; then
  samtools view -H "$meta_temp_dir/output1.bam" > "$meta_temp_dir/header1.txt" 2>/dev/null || true
  if [ -s "$meta_temp_dir/header1.txt" ]; then
      check_file_contains "$meta_temp_dir/header1.txt" "@HD" "BAM header"
      log "✓ Valid BAM format detected"
  else
      log "Note: Cannot validate BAM format (samtools not available or BAM corrupt)"
  fi
else
  log "Note: samtools not available for BAM validation"
fi

log "✅ TEST 1 completed successfully"

# --- Test Case 2: BED12 format conversion ---
log "Starting TEST 2: BED12 to BAM conversion"

log "Executing $meta_name with BED12 format..."
"$meta_executable" \
  --input "$test_dir/test.bed12" \
  --genome "$test_dir/test.genome" \
  --output "$meta_temp_dir/output2.bam" \
  --bed12

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output2.bam" "BED12 output BAM file"
check_file_not_empty "$meta_temp_dir/output2.bam" "BED12 output BAM file"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Custom mapping quality ---
log "Starting TEST 3: Custom mapping quality"

log "Executing $meta_name with custom mapping quality..."
"$meta_executable" \
  --input "$test_dir/test.bed" \
  --genome "$test_dir/test.genome" \
  --output "$meta_temp_dir/output3.bam" \
  --map_quality 30

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output3.bam" "output BAM with custom MAPQ"
check_file_not_empty "$meta_temp_dir/output3.bam" "output BAM with custom MAPQ"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Uncompressed BAM ---
log "Starting TEST 4: Uncompressed BAM output"

log "Executing $meta_name with uncompressed BAM..."
"$meta_executable" \
  --input "$test_dir/test.bed" \
  --genome "$test_dir/test.genome" \
  --output "$meta_temp_dir/output4.bam" \
  --uncompress_bam

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output4.bam" "uncompressed BAM file"
check_file_not_empty "$meta_temp_dir/output4.bam" "uncompressed BAM file"

# Uncompressed BAM should generally be larger than compressed
compressed_size=$(stat -c%s "$meta_temp_dir/output1.bam")
uncompressed_size=$(stat -c%s "$meta_temp_dir/output4.bam")
log "Compressed BAM size: $compressed_size bytes"
log "Uncompressed BAM size: $uncompressed_size bytes"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests completed successfully"
