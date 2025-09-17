#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_shift"

####################################################################################################

log "Creating test data..."

# Create test BED file with intervals
cat > "$meta_temp_dir/input.bed" << 'EOF'
chr1	100	200	feature1	100	+
chr1	300	400	feature2	200	-
chr2	150	250	feature3	150	+
chr2	350	450	feature4	300	-
chr3	500	600	feature5	400	+
EOF

# Create genome file defining chromosome sizes
cat > "$meta_temp_dir/genome.txt" << 'EOF'
chr1	1000
chr2	1000
chr3	1000
EOF

check_file_exists "$meta_temp_dir/input.bed" "input BED file"
check_file_exists "$meta_temp_dir/genome.txt" "genome file"

####################################################################################################

log "TEST 1: Basic shift functionality with positive shift"

"$meta_executable" \
  --input "$meta_temp_dir/input.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --shift 50 \
  --output "$meta_temp_dir/shift_positive.bed"

check_file_exists "$meta_temp_dir/shift_positive.bed" "positive shift output"
check_file_not_empty "$meta_temp_dir/shift_positive.bed" "positive shift output"

# Check that coordinates were shifted correctly
# First feature should be shifted from 100-200 to 150-250
if grep -q "chr1	150	250	feature1" "$meta_temp_dir/shift_positive.bed"; then
  log "✓ Positive shift correctly applied to first feature"
else
  log "✗ Positive shift not correctly applied"
  cat "$meta_temp_dir/shift_positive.bed"
  exit 1
fi

####################################################################################################

log "TEST 2: Basic shift functionality with negative shift"

"$meta_executable" \
  --input "$meta_temp_dir/input.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --shift -25 \
  --output "$meta_temp_dir/shift_negative.bed"

check_file_exists "$meta_temp_dir/shift_negative.bed" "negative shift output"

# Check that coordinates were shifted correctly
# First feature should be shifted from 100-200 to 75-175
if grep -q "chr1	75	175	feature1" "$meta_temp_dir/shift_negative.bed"; then
  log "✓ Negative shift correctly applied"
else
  log "✗ Negative shift not correctly applied"
  cat "$meta_temp_dir/shift_negative.bed"
  exit 1
fi

####################################################################################################

log "TEST 3: Strand-specific shifting"

"$meta_executable" \
  --input "$meta_temp_dir/input.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --plus_shift 30 \
  --minus_shift -20 \
  --output "$meta_temp_dir/shift_strand.bed"

check_file_exists "$meta_temp_dir/shift_strand.bed" "strand-specific shift output"

# Check plus strand feature (feature1: + strand)
if grep -q "chr1	130	230	feature1.*+" "$meta_temp_dir/shift_strand.bed"; then
  log "✓ Plus strand shift correctly applied"
else
  log "✗ Plus strand shift not correctly applied"
  cat "$meta_temp_dir/shift_strand.bed"
  exit 1
fi

# Check minus strand feature (feature2: - strand)
if grep -q "chr1	280	380	feature2.*-" "$meta_temp_dir/shift_strand.bed"; then
  log "✓ Minus strand shift correctly applied"
else
  log "✗ Minus strand shift not correctly applied"
  cat "$meta_temp_dir/shift_strand.bed"
  exit 1
fi

####################################################################################################

log "TEST 4: Boundary handling - shift below 0"

# Create interval close to start of chromosome
cat > "$meta_temp_dir/boundary.bed" << 'EOF'
chr1	10	50	near_start	100	+
EOF

"$meta_executable" \
  --input "$meta_temp_dir/boundary.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --shift -20 \
  --output "$meta_temp_dir/boundary_test.bed"

check_file_exists "$meta_temp_dir/boundary_test.bed" "boundary test output"

# Start should be clamped to 0, end should be 30 (original 50 - 20)
if grep -q "chr1	0	30" "$meta_temp_dir/boundary_test.bed"; then
  log "✓ Boundary clamping to 0 works correctly"
else
  log "✗ Boundary clamping failed"
  cat "$meta_temp_dir/boundary_test.bed"
  exit 1
fi

####################################################################################################

log "TEST 5: Boundary handling - shift beyond chromosome length"

# Create interval near end of chromosome
cat > "$meta_temp_dir/boundary_end.bed" << 'EOF'
chr1	950	990	near_end	100	+
EOF

"$meta_executable" \
  --input "$meta_temp_dir/boundary_end.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --shift 50 \
  --output "$meta_temp_dir/boundary_end_test.bed"

check_file_exists "$meta_temp_dir/boundary_end_test.bed" "boundary end test output"

# End should be clamped to chromosome length (1000), start adjusted to maintain valid interval
if grep -q "chr1	999	1000" "$meta_temp_dir/boundary_end_test.bed"; then
  log "✓ Boundary clamping to chromosome length works correctly"
else
  log "✗ End boundary clamping failed"
  cat "$meta_temp_dir/boundary_end_test.bed"
  exit 1
fi

####################################################################################################

log "TEST 6: Percentage-based shifting"

# Create test interval of known length (100bp)
cat > "$meta_temp_dir/percentage.bed" << 'EOF'
chr1	100	200	test_feature	100	+
EOF

"$meta_executable" \
  --input "$meta_temp_dir/percentage.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --shift 0.5 \
  --percentage \
  --output "$meta_temp_dir/percentage_test.bed"

check_file_exists "$meta_temp_dir/percentage_test.bed" "percentage test output"

# 50% of 100bp = 50bp shift, so 100-200 becomes 150-250
if grep -q "chr1	150	250" "$meta_temp_dir/percentage_test.bed"; then
  log "✓ Percentage-based shifting works correctly"
else
  log "✗ Percentage-based shifting failed"
  cat "$meta_temp_dir/percentage_test.bed"
  exit 1
fi

####################################################################################################

log "TEST 7: Header preservation"

# Create input with header
cat > "$meta_temp_dir/with_header.bed" << 'EOF'
# BED file header
# Track information
chr1	100	200	feature1	100	+
chr1	300	400	feature2	200	-
EOF

"$meta_executable" \
  --input "$meta_temp_dir/with_header.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --shift 10 \
  --header \
  --output "$meta_temp_dir/header_test.bed"

check_file_exists "$meta_temp_dir/header_test.bed" "header test output"

# Check that header lines are preserved
if grep -q "# BED file header" "$meta_temp_dir/header_test.bed"; then
  log "✓ Header preservation works correctly"
else
  log "✗ Header preservation failed"
  cat "$meta_temp_dir/header_test.bed"
  exit 1
fi

####################################################################################################

log "TEST 8: Parameter validation"

# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" \
  --genome "$meta_temp_dir/genome.txt" \
  --shift 10 \
  --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "✗ Should have failed without --input parameter"
  exit 1
else
  log "✓ Correctly requires --input parameter"
fi

if "$meta_executable" \
  --input "$meta_temp_dir/input.bed" \
  --shift 10 \
  --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "✗ Should have failed without --genome parameter"
  exit 1
else
  log "✓ Correctly requires --genome parameter"
fi

if "$meta_executable" \
  --input "$meta_temp_dir/input.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "✗ Should have failed without shift parameters"
  exit 1
else
  log "✓ Correctly requires shift parameters"
fi

####################################################################################################

log "TEST 9: Parameter combination validation"

# Test invalid parameter combinations
log "Testing parameter combination validation"

if "$meta_executable" \
  --input "$meta_temp_dir/input.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --shift 10 \
  --plus_shift 5 \
  --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "✗ Should have failed with conflicting shift parameters"
  exit 1
else
  log "✓ Correctly rejects conflicting shift parameters"
fi

if "$meta_executable" \
  --input "$meta_temp_dir/input.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --plus_shift 5 \
  --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "✗ Should have failed with incomplete strand-specific parameters"
  exit 1
else
  log "✓ Correctly requires both plus and minus shift parameters"
fi

####################################################################################################

log "TEST 10: File validation"

# Test with non-existent files
if "$meta_executable" \
  --input "/nonexistent/file.bed" \
  --genome "$meta_temp_dir/genome.txt" \
  --shift 10 \
  --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "Should have failed with non-existent input file"
else
  log "✓ Properly handles non-existent input files"
fi

####################################################################################################

log "All tests completed successfully!"
