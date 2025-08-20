#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_slop"

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

log "TEST 1: Symmetric extension with --both parameter"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both 50 \
    --output "$meta_temp_dir/slop_both.bed"

check_file_exists "$meta_temp_dir/slop_both.bed" "symmetric extension output"
check_file_not_empty "$meta_temp_dir/slop_both.bed" "symmetric extension output"

# Check that coordinates were extended correctly
# First feature should be extended from 100-200 to 50-250
if grep -q "chr1	50	250	feature1" "$meta_temp_dir/slop_both.bed"; then
    log "✓ Symmetric extension correctly applied to first feature"
else
    log "✗ Symmetric extension not correctly applied"
    cat "$meta_temp_dir/slop_both.bed"
    exit 1
fi

####################################################################################################

log "TEST 2: Asymmetric extension with --left and --right parameters"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --left 30 \
    --right 20 \
    --output "$meta_temp_dir/slop_asymmetric.bed"

check_file_exists "$meta_temp_dir/slop_asymmetric.bed" "asymmetric extension output"

# Check that coordinates were extended correctly
# First feature should be extended from 100-200 to 70-220 (start-30, end+20)
if grep -q "chr1	70	220	feature1" "$meta_temp_dir/slop_asymmetric.bed"; then
    log "✓ Asymmetric extension correctly applied"
else
    log "✗ Asymmetric extension not correctly applied"
    cat "$meta_temp_dir/slop_asymmetric.bed"
    exit 1
fi

####################################################################################################

log "TEST 3: Strand-aware extension"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --left 40 \
    --right 10 \
    --strand_aware \
    --output "$meta_temp_dir/slop_strand.bed"

check_file_exists "$meta_temp_dir/slop_strand.bed" "strand-aware extension output"

# For strand-aware mode:
# + strand: left extends start (upstream), right extends end (downstream)
# - strand: left extends end (downstream), right extends start (upstream)

# Check plus strand feature (feature1: + strand) - should be like normal: start-40, end+10
if grep -q "chr1	60	210	feature1.*+" "$meta_temp_dir/slop_strand.bed"; then
    log "✓ Plus strand extension correctly applied"
else
    log "✗ Plus strand extension not correctly applied"
    cat "$meta_temp_dir/slop_strand.bed"
    exit 1
fi

# Check minus strand feature (feature2: - strand) - should be: start-10, end+40
if grep -q "chr1	290	440	feature2.*-" "$meta_temp_dir/slop_strand.bed"; then
    log "✓ Minus strand extension correctly applied"
else
    log "✗ Minus strand extension not correctly applied"
    cat "$meta_temp_dir/slop_strand.bed"
    exit 1
fi

####################################################################################################

log "TEST 4: Boundary handling - extension below 0"

# Create interval close to start of chromosome
cat > "$meta_temp_dir/boundary.bed" << 'EOF'
chr1	10	50	near_start	100	+
EOF

"$meta_executable" \
    --input "$meta_temp_dir/boundary.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both 20 \
    --output "$meta_temp_dir/boundary_test.bed"

check_file_exists "$meta_temp_dir/boundary_test.bed" "boundary test output"

# Start should be clamped to 0, end should be 70 (original 50 + 20)
if grep -q "chr1	0	70" "$meta_temp_dir/boundary_test.bed"; then
    log "✓ Boundary clamping to 0 works correctly"
else
    log "✗ Boundary clamping failed"
    cat "$meta_temp_dir/boundary_test.bed"
    exit 1
fi

####################################################################################################

log "TEST 5: Boundary handling - extension beyond chromosome length"

# Create interval near end of chromosome
cat > "$meta_temp_dir/boundary_end.bed" << 'EOF'
chr1	950	980	near_end	100	+
EOF

"$meta_executable" \
    --input "$meta_temp_dir/boundary_end.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both 50 \
    --output "$meta_temp_dir/boundary_end_test.bed"

check_file_exists "$meta_temp_dir/boundary_end_test.bed" "boundary end test output"

# Start should be 900 (950-50), end should be clamped to chromosome length (1000)
if grep -q "chr1	900	1000" "$meta_temp_dir/boundary_end_test.bed"; then
    log "✓ Boundary clamping to chromosome length works correctly"
else
    log "✗ End boundary clamping failed"
    cat "$meta_temp_dir/boundary_end_test.bed"
    exit 1
fi

####################################################################################################

log "TEST 6: Percentage-based extension"

# Create test interval of known length (100bp)
cat > "$meta_temp_dir/percentage.bed" << 'EOF'
chr1	100	200	test_feature	100	+
EOF

"$meta_executable" \
    --input "$meta_temp_dir/percentage.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both 0.5 \
    --percentage \
    --output "$meta_temp_dir/percentage_test.bed"

check_file_exists "$meta_temp_dir/percentage_test.bed" "percentage test output"

# 50% of 100bp = 50bp extension in each direction, so 100-200 becomes 50-250
if grep -q "chr1	50	250" "$meta_temp_dir/percentage_test.bed"; then
    log "✓ Percentage-based extension works correctly"
else
    log "✗ Percentage-based extension failed"
    cat "$meta_temp_dir/percentage_test.bed"
    exit 1
fi

####################################################################################################

log "TEST 7: Asymmetric percentage-based extension"

"$meta_executable" \
    --input "$meta_temp_dir/percentage.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --left 0.3 \
    --right 0.2 \
    --percentage \
    --output "$meta_temp_dir/percentage_asymmetric.bed"

check_file_exists "$meta_temp_dir/percentage_asymmetric.bed" "asymmetric percentage test output"

# 30% of 100bp = 30bp left, 20% of 100bp = 20bp right, so 100-200 becomes 70-220
if grep -q "chr1	70	220" "$meta_temp_dir/percentage_asymmetric.bed"; then
    log "✓ Asymmetric percentage-based extension works correctly"
else
    log "✗ Asymmetric percentage-based extension failed"
    cat "$meta_temp_dir/percentage_asymmetric.bed"
    exit 1
fi

####################################################################################################

log "TEST 8: Header preservation"

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
    --both 10 \
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

log "TEST 9: Parameter validation"

# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" \
    --genome "$meta_temp_dir/genome.txt" \
    --both 10 \
    --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "✗ Should have failed without --input parameter"
    exit 1
else
    log "✓ Correctly requires --input parameter"
fi

if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --both 10 \
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
    log "✗ Should have failed without extension parameters"
    exit 1
else
    log "✓ Correctly requires extension parameters"
fi

####################################################################################################

log "TEST 10: Parameter combination validation"

# Test invalid parameter combinations
log "Testing parameter combination validation"

if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both 10 \
    --left 5 \
    --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "✗ Should have failed with conflicting extension parameters"
    exit 1
else
    log "✓ Correctly rejects conflicting extension parameters"
fi

if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --left 5 \
    --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "✗ Should have failed with incomplete left/right parameters"
    exit 1
else
    log "✓ Correctly requires both left and right parameters together"
fi

####################################################################################################

log "TEST 11: File validation"

# Test with non-existent files
if "$meta_executable" \
    --input "/nonexistent/file.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both 10 \
    --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "Should have failed with non-existent input file"
else
    log "✓ Properly handles non-existent input files"
fi

####################################################################################################

log "TEST 12: Zero extension test"

# Test with zero extension (should preserve original coordinates)
"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both 0 \
    --output "$meta_temp_dir/zero_extension.bed"

check_file_exists "$meta_temp_dir/zero_extension.bed" "zero extension output"

# Check that first feature coordinates are unchanged
if grep -q "chr1	100	200	feature1" "$meta_temp_dir/zero_extension.bed"; then
    log "✓ Zero extension preserves original coordinates"
else
    log "✗ Zero extension failed"
    cat "$meta_temp_dir/zero_extension.bed"
    exit 1
fi

####################################################################################################

log "All tests completed successfully!"
