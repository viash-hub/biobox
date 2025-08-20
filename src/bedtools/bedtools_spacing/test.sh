#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_spacing"

####################################################################################################

log "Creating test data..."

# Create test BED file with intervals demonstrating different spacing scenarios
# Note: bedtools spacing requires sorted input by chromosome and start coordinate
cat > "$meta_temp_dir/input.bed" << 'EOF'
chr1	0	10	feature1	100	+
chr1	10	20	feature2	200	-
chr1	19	30	feature3	150	+
chr1	35	45	feature4	300	-
chr1	100	200	feature5	400	+
chr2	50	60	feature6	100	+
chr2	70	80	feature7	200	-
chr2	85	95	feature8	150	+
EOF

check_file_exists "$meta_temp_dir/input.bed" "input BED file"

####################################################################################################

log "TEST 1: Basic spacing analysis"

"$meta_executable" \
  --input "$meta_temp_dir/input.bed" \
  --output "$meta_temp_dir/spacing_basic.bed"

check_file_exists "$meta_temp_dir/spacing_basic.bed" "basic spacing output"
check_file_not_empty "$meta_temp_dir/spacing_basic.bed" "basic spacing output"

# Check that output has the same number of lines as input
input_lines=$(wc -l < "$meta_temp_dir/input.bed")
output_lines=$(wc -l < "$meta_temp_dir/spacing_basic.bed")

if [ "$input_lines" -eq "$output_lines" ]; then
  log "✓ Same number of intervals preserved ($output_lines)"
else
  log "✗ Number of intervals changed: input=$input_lines, output=$output_lines"
  exit 1
fi

# Verify that spacing column was added (should have one more column than input)
input_cols=$(head -1 "$meta_temp_dir/input.bed" | awk '{print NF}')
output_cols=$(head -1 "$meta_temp_dir/spacing_basic.bed" | awk '{print NF}')

if [ "$output_cols" -eq $((input_cols + 1)) ]; then
  log "✓ Spacing column added correctly"
else
  log "✗ Expected $((input_cols + 1)) columns, got $output_cols"
  exit 1
fi

####################################################################################################

log "TEST 2: Verify specific spacing calculations"

# Check specific spacing patterns from our test data:
# chr1 0-10: first interval, should have "."
# chr1 10-20: adjacent (0 gap) to previous
# chr1 19-30: overlapping (-1) with previous 
# chr1 35-45: gap of 5 bp from previous
# chr1 100-200: gap of 55 bp from previous

# First interval on chr1 should have "." (NULL distance)
if grep -q "chr1	0	10	feature1	100	+	\." "$meta_temp_dir/spacing_basic.bed"; then
  log "✓ First interval correctly marked with NULL distance (.)"
else
  log "✗ First interval should have NULL distance"
  cat "$meta_temp_dir/spacing_basic.bed"
  exit 1
fi

# Second interval should be adjacent (0 gap)
if grep -q "chr1	10	20	feature2	200	-	0" "$meta_temp_dir/spacing_basic.bed"; then
  log "✓ Adjacent intervals correctly show 0 gap"
else
  log "✗ Adjacent intervals should show 0 gap"
  cat "$meta_temp_dir/spacing_basic.bed"
  exit 1
fi

# Third interval overlaps with second (-1)
if grep -q "chr1	19	30	feature3	150	+	-1" "$meta_temp_dir/spacing_basic.bed"; then
  log "✓ Overlapping intervals correctly show -1"
else
  log "✗ Overlapping intervals should show -1"
  cat "$meta_temp_dir/spacing_basic.bed"
  exit 1
fi

# Fourth interval has gap of 5 bp (35 - 30 = 5)
if grep -q "chr1	35	45	feature4	300	-	5" "$meta_temp_dir/spacing_basic.bed"; then
  log "✓ Gap spacing correctly calculated (5 bp)"
else
  log "✗ Gap spacing calculation incorrect"
  cat "$meta_temp_dir/spacing_basic.bed"
  exit 1
fi

# Fifth interval has gap of 55 bp (100 - 45 = 55)
if grep -q "chr1	100	200	feature5	400	+	55" "$meta_temp_dir/spacing_basic.bed"; then
  log "✓ Larger gap correctly calculated (55 bp)"
else
  log "✗ Larger gap calculation incorrect"
  cat "$meta_temp_dir/spacing_basic.bed"
  exit 1
fi

####################################################################################################

log "TEST 3: Multiple chromosomes handling"

# Check that first interval on chr2 gets NULL distance
if grep -q "chr2	50	60	feature6	100	+	\." "$meta_temp_dir/spacing_basic.bed"; then
  log "✓ First interval on chr2 correctly marked with NULL distance"
else
  log "✗ First interval on new chromosome should have NULL distance"
  cat "$meta_temp_dir/spacing_basic.bed"
  exit 1
fi

# Check gap calculation within chr2 (70-60=10)
if grep -q "chr2	70	80	feature7	200	-	10" "$meta_temp_dir/spacing_basic.bed"; then
  log "✓ Gap within chr2 correctly calculated (10 bp)"
else
  log "✗ Gap calculation within chr2 incorrect"
  cat "$meta_temp_dir/spacing_basic.bed"
  exit 1
fi

# Check adjacent intervals on chr2 (85-80=5)
if grep -q "chr2	85	95	feature8	150	+	5" "$meta_temp_dir/spacing_basic.bed"; then
  log "✓ Gap on chr2 correctly calculated (5 bp)"
else
  log "✗ Second gap calculation within chr2 incorrect"
  cat "$meta_temp_dir/spacing_basic.bed"
  exit 1
fi

####################################################################################################

log "TEST 4: Header preservation"

# Create input with header
cat > "$meta_temp_dir/with_header.bed" << 'EOF'
# BED file header
# Track information
chr1	100	200	feature1	100	+
chr1	250	350	feature2	200	-
EOF

"$meta_executable" \
  --input "$meta_temp_dir/with_header.bed" \
  --output "$meta_temp_dir/header_test.bed" \
  --include_header

check_file_exists "$meta_temp_dir/header_test.bed" "header test output"

# Check that header lines are preserved
if grep -q "# BED file header" "$meta_temp_dir/header_test.bed"; then
  log "✓ Header preservation works correctly"
else
  log "✗ Header preservation failed"
  cat "$meta_temp_dir/header_test.bed"
  exit 1
fi

# Verify spacing calculation still works with header
if grep -q "chr1	100	200	feature1	100	+	\." "$meta_temp_dir/header_test.bed"; then
  log "✓ Spacing calculation works with header present"
else
  log "✗ Spacing calculation failed with header"
  cat "$meta_temp_dir/header_test.bed"
  exit 1
fi

####################################################################################################

log "TEST 5: Adjacent intervals (0 gap) test"

# Create test with perfectly adjacent intervals
cat > "$meta_temp_dir/adjacent.bed" << 'EOF'
chr1	100	200	interval1	100	+
chr1	200	300	interval2	200	+
chr1	300	400	interval3	300	+
EOF

"$meta_executable" \
  --input "$meta_temp_dir/adjacent.bed" \
  --output "$meta_temp_dir/adjacent_test.bed"

check_file_exists "$meta_temp_dir/adjacent_test.bed" "adjacent intervals test"

# All intervals after the first should have 0 gap
if grep -q "chr1	200	300	interval2	200	+	0" "$meta_temp_dir/adjacent_test.bed" && \
   grep -q "chr1	300	400	interval3	300	+	0" "$meta_temp_dir/adjacent_test.bed"; then
  log "✓ Adjacent intervals correctly show 0 gaps"
else
  log "✗ Adjacent intervals spacing calculation failed"
  cat "$meta_temp_dir/adjacent_test.bed"
  exit 1
fi

####################################################################################################

log "TEST 6: Overlapping intervals test"

# Create test with overlapping intervals
cat > "$meta_temp_dir/overlapping.bed" << 'EOF'
chr1	100	200	interval1	100	+
chr1	150	250	interval2	200	+
chr1	180	280	interval3	300	+
EOF

"$meta_executable" \
  --input "$meta_temp_dir/overlapping.bed" \
  --output "$meta_temp_dir/overlapping_test.bed"

check_file_exists "$meta_temp_dir/overlapping_test.bed" "overlapping intervals test"

# All overlapping intervals should have -1
if grep -q "chr1	150	250	interval2	200	+	-1" "$meta_temp_dir/overlapping_test.bed" && \
   grep -q "chr1	180	280	interval3	300	+	-1" "$meta_temp_dir/overlapping_test.bed"; then
  log "✓ Overlapping intervals correctly show -1 gaps"
else
  log "✗ Overlapping intervals spacing calculation failed"
  cat "$meta_temp_dir/overlapping_test.bed"
  exit 1
fi

####################################################################################################

log "TEST 7: Parameter validation"

# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" \
  --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "✗ Should have failed without --input parameter"
  exit 1
else
  log "✓ Correctly requires --input parameter"
fi

if "$meta_executable" \
  --input "$meta_temp_dir/input.bed" 2>/dev/null; then
  log "✗ Should have failed without --output parameter"
  exit 1
else
  log "✓ Correctly requires --output parameter"
fi

####################################################################################################

log "TEST 8: File validation"

# Test with non-existent files
if "$meta_executable" \
  --input "/nonexistent/file.bed" \
  --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "Should have failed with non-existent input file"
else
  log "✓ Properly handles non-existent input files"
fi

####################################################################################################

log "TEST 9: Empty input handling"

# Create empty input file
touch "$meta_temp_dir/empty.bed"

"$meta_executable" \
  --input "$meta_temp_dir/empty.bed" \
  --output "$meta_temp_dir/empty_output.bed"

check_file_exists "$meta_temp_dir/empty_output.bed" "empty input test output"

# Output should also be empty
if [ ! -s "$meta_temp_dir/empty_output.bed" ]; then
  log "✓ Empty input produces empty output"
else
  log "✗ Empty input handling failed"
  cat "$meta_temp_dir/empty_output.bed"
  exit 1
fi

####################################################################################################

log "TEST 10: Single interval test"

# Create test with single interval
cat > "$meta_temp_dir/single.bed" << 'EOF'
chr1	100	200	single_interval	100	+
EOF

"$meta_executable" \
  --input "$meta_temp_dir/single.bed" \
  --output "$meta_temp_dir/single_test.bed"

check_file_exists "$meta_temp_dir/single_test.bed" "single interval test"

# Single interval should have NULL distance
if grep -q "chr1	100	200	single_interval	100	+	\." "$meta_temp_dir/single_test.bed"; then
  log "✓ Single interval correctly shows NULL distance"
else
  log "✗ Single interval spacing failed"
  cat "$meta_temp_dir/single_test.bed"
  exit 1
fi

####################################################################################################

log "TEST 11: Performance options test"

# Test with buffer options (basic functionality test)
"$meta_executable" \
  --input "$meta_temp_dir/input.bed" \
  --output "$meta_temp_dir/nobuf_test.bed" \
  --no_buffer

check_file_exists "$meta_temp_dir/nobuf_test.bed" "no buffer test output"

# Should produce same results as buffered version
if diff -q "$meta_temp_dir/spacing_basic.bed" "$meta_temp_dir/nobuf_test.bed" >/dev/null; then
  log "✓ No buffer option produces identical results"
else
  log "✗ No buffer option changed results"
  exit 1
fi

####################################################################################################

log "All tests completed successfully!"
