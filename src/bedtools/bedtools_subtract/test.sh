#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_subtract"

####################################################################################################

log "Creating test data..."

# Create test BED file A (target intervals)
cat > "$meta_temp_dir/input_a.bed" << 'EOF'
chr1	100	200	feature1	100	+
chr1	300	400	feature2	200	+
chr1	500	600	feature3	300	+
chr1	700	800	feature4	400	-
chr2	100	300	feature5	500	+
chr2	400	500	feature6	600	-
EOF

# Create test BED file B (intervals to subtract)
cat > "$meta_temp_dir/input_b.bed" << 'EOF'
chr1	150	175	repeat1	50	.
chr1	350	450	repeat2	75	.
chr1	550	650	repeat3	100	.
chr2	120	180	repeat4	60	.
chr2	420	480	repeat5	80	.
EOF

check_file_exists "$meta_temp_dir/input_a.bed" "input file A"
check_file_exists "$meta_temp_dir/input_b.bed" "input file B"

####################################################################################################

log "TEST 1: Basic subtraction (partial overlap removal)"

"$meta_executable" \
    --input_a "$meta_temp_dir/input_a.bed" \
    --input_b "$meta_temp_dir/input_b.bed" \
    --output "$meta_temp_dir/basic_subtract.bed"

check_file_exists "$meta_temp_dir/basic_subtract.bed" "basic subtraction output"
check_file_not_empty "$meta_temp_dir/basic_subtract.bed" "basic subtraction output"

# Check that overlapping portions were removed and non-overlapping parts preserved
# chr1 100-200 overlaps with 150-175, should create two intervals: 100-150 and 175-200
if grep -q "chr1	100	150	feature1	100	+" "$meta_temp_dir/basic_subtract.bed" && \
   grep -q "chr1	175	200	feature1	100	+" "$meta_temp_dir/basic_subtract.bed"; then
    log "✓ Partial overlap correctly creates split intervals"
else
    log "✗ Partial overlap handling failed"
    cat "$meta_temp_dir/basic_subtract.bed"
    exit 1
fi

# Check that feature2 (300-400) overlapping with repeat2 (350-450) creates interval 300-350
if grep -q "chr1	300	350	feature2	200	+" "$meta_temp_dir/basic_subtract.bed"; then
    log "✓ Partial overlap at end correctly handled"
else
    log "✗ End overlap handling failed"
    cat "$meta_temp_dir/basic_subtract.bed"
    exit 1
fi

####################################################################################################

log "TEST 2: Complete overlap removal"

# Create test where B completely covers some A intervals
cat > "$meta_temp_dir/input_a2.bed" << 'EOF'
chr1	100	200	small1	100	+
chr1	300	400	small2	200	+
chr1	500	600	small3	300	+
EOF

cat > "$meta_temp_dir/input_b2.bed" << 'EOF'
chr1	50	250	big1	100	.
chr1	550	580	partial1	50	.
EOF

"$meta_executable" \
    --input_a "$meta_temp_dir/input_a2.bed" \
    --input_b "$meta_temp_dir/input_b2.bed" \
    --output "$meta_temp_dir/complete_overlap.bed"

check_file_exists "$meta_temp_dir/complete_overlap.bed" "complete overlap output"

# First interval should be completely removed (covered by big1)
# Second interval should remain (no overlap)
# Third interval should have partial removal (550-580 removed from 500-600)
if ! grep -q "chr1	100	200" "$meta_temp_dir/complete_overlap.bed" && \
   grep -q "chr1	300	400	small2	200	+" "$meta_temp_dir/complete_overlap.bed" && \
   grep -q "chr1	500	550	small3	300	+" "$meta_temp_dir/complete_overlap.bed" && \
   grep -q "chr1	580	600	small3	300	+" "$meta_temp_dir/complete_overlap.bed"; then
    log "✓ Complete and partial overlap handling works correctly"
else
    log "✗ Complete/partial overlap handling failed"
    cat "$meta_temp_dir/complete_overlap.bed"
    exit 1
fi

####################################################################################################

log "TEST 3: Remove entire feature (-A option)"

"$meta_executable" \
    --input_a "$meta_temp_dir/input_a.bed" \
    --input_b "$meta_temp_dir/input_b.bed" \
    --output "$meta_temp_dir/remove_entire.bed" \
    --remove_entire

check_file_exists "$meta_temp_dir/remove_entire.bed" "remove entire output"

# With -A, any interval in A that overlaps with B should be completely removed
# Only intervals with no overlap should remain
if ! grep -q "chr1	100	200" "$meta_temp_dir/remove_entire.bed" && \
   ! grep -q "chr1	300	400" "$meta_temp_dir/remove_entire.bed" && \
   ! grep -q "chr1	500	600" "$meta_temp_dir/remove_entire.bed" && \
   ! grep -q "chr2	100	300" "$meta_temp_dir/remove_entire.bed" && \
   ! grep -q "chr2	400	500" "$meta_temp_dir/remove_entire.bed" && \
   grep -q "chr1	700	800	feature4	400	-" "$meta_temp_dir/remove_entire.bed"; then
    log "✓ Remove entire feature (-A) works correctly"
else
    log "✗ Remove entire feature (-A) failed"
    cat "$meta_temp_dir/remove_entire.bed"
    exit 1
fi

####################################################################################################

log "TEST 4: Minimum overlap fraction (-f option)"

"$meta_executable" \
    --input_a "$meta_temp_dir/input_a.bed" \
    --input_b "$meta_temp_dir/input_b.bed" \
    --output "$meta_temp_dir/min_overlap.bed" \
    --min_overlap_a 0.5

check_file_exists "$meta_temp_dir/min_overlap.bed" "minimum overlap output"

# With -f 0.5, B must overlap at least 50% of A for subtraction to occur
# chr1 100-200 (length 100): repeat1 150-175 (length 25) overlaps 25bp = 25% < 50%, no subtraction
# chr1 300-400 (length 100): repeat2 350-450 overlaps 50bp = 50% >= 50%, subtraction occurs
if grep -q "chr1	100	200	feature1	100	+" "$meta_temp_dir/min_overlap.bed" && \
   grep -q "chr1	300	350	feature2	200	+" "$meta_temp_dir/min_overlap.bed"; then
    log "✓ Minimum overlap fraction (-f) works correctly"
else
    log "✗ Minimum overlap fraction (-f) failed"
    cat "$meta_temp_dir/min_overlap.bed"
    exit 1
fi

####################################################################################################

log "TEST 5: Strand-specific subtraction (-s option)"

# Create strand-specific test data
cat > "$meta_temp_dir/strand_a.bed" << 'EOF'
chr1	100	200	plus1	100	+
chr1	300	400	minus1	200	-
chr1	500	600	plus2	300	+
EOF

cat > "$meta_temp_dir/strand_b.bed" << 'EOF'
chr1	150	175	plus_repeat	50	+
chr1	350	375	minus_repeat	75	-
chr1	550	575	minus_repeat2	100	-
EOF

"$meta_executable" \
    --input_a "$meta_temp_dir/strand_a.bed" \
    --input_b "$meta_temp_dir/strand_b.bed" \
    --output "$meta_temp_dir/strand_specific.bed" \
    --same_strand

check_file_exists "$meta_temp_dir/strand_specific.bed" "strand-specific output"

# With -s, only same-strand overlaps cause subtraction
# plus1 (+) overlaps with plus_repeat (+) -> subtraction occurs
# minus1 (-) overlaps with minus_repeat (-) -> subtraction occurs  
# plus2 (+) overlaps with minus_repeat2 (-) -> no subtraction (different strands)
if grep -q "chr1	100	150	plus1	100	+" "$meta_temp_dir/strand_specific.bed" && \
   grep -q "chr1	175	200	plus1	100	+" "$meta_temp_dir/strand_specific.bed" && \
   grep -q "chr1	300	350	minus1	200	-" "$meta_temp_dir/strand_specific.bed" && \
   grep -q "chr1	375	400	minus1	200	-" "$meta_temp_dir/strand_specific.bed" && \
   grep -q "chr1	500	600	plus2	300	+" "$meta_temp_dir/strand_specific.bed"; then
    log "✓ Strand-specific subtraction (-s) works correctly"
else
    log "✗ Strand-specific subtraction (-s) failed"
    cat "$meta_temp_dir/strand_specific.bed"
    exit 1
fi

####################################################################################################

log "TEST 6: Write overlap information (-wo option)"

"$meta_executable" \
    --input_a "$meta_temp_dir/input_a.bed" \
    --input_b "$meta_temp_dir/input_b.bed" \
    --output "$meta_temp_dir/overlap_info.bed" \
    --write_overlap_counts

check_file_exists "$meta_temp_dir/overlap_info.bed" "overlap info output"
check_file_not_empty "$meta_temp_dir/overlap_info.bed" "overlap info output"

# With -wo, output includes A interval, B interval, and overlap count
# Should have more columns than basic subtraction
input_cols=$(head -1 "$meta_temp_dir/input_a.bed" | awk '{print NF}')
output_cols=$(head -1 "$meta_temp_dir/overlap_info.bed" | awk '{print NF}')

if [ "$output_cols" -gt "$input_cols" ]; then
    log "✓ Overlap information (-wo) adds additional columns"
else
    log "✗ Overlap information (-wo) format incorrect"
    cat "$meta_temp_dir/overlap_info.bed"
    exit 1
fi

####################################################################################################

log "TEST 7: Header preservation"

# Create input with header
cat > "$meta_temp_dir/with_header_a.bed" << 'EOF'
# BED file header A
# Track information
chr1	100	200	feature1	100	+
chr1	300	400	feature2	200	+
EOF

cat > "$meta_temp_dir/with_header_b.bed" << 'EOF'
chr1	150	175	repeat1	50	.
EOF

"$meta_executable" \
    --input_a "$meta_temp_dir/with_header_a.bed" \
    --input_b "$meta_temp_dir/with_header_b.bed" \
    --output "$meta_temp_dir/header_test.bed" \
    --include_header

check_file_exists "$meta_temp_dir/header_test.bed" "header test output"

# Check that header lines are preserved
if grep -q "# BED file header A" "$meta_temp_dir/header_test.bed"; then
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
    --input_b "$meta_temp_dir/input_b.bed" \
    --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "✗ Should have failed without --input_a parameter"
    exit 1
else
    log "✓ Correctly requires --input_a parameter"
fi

if "$meta_executable" \
    --input_a "$meta_temp_dir/input_a.bed" \
    --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "✗ Should have failed without --input_b parameter"
    exit 1
else
    log "✓ Correctly requires --input_b parameter"
fi

if "$meta_executable" \
    --input_a "$meta_temp_dir/input_a.bed" \
    --input_b "$meta_temp_dir/input_b.bed" 2>/dev/null; then
    log "✗ Should have failed without --output parameter"
    exit 1
else
    log "✓ Correctly requires --output parameter"
fi

####################################################################################################

log "TEST 9: File validation"

# Test with non-existent files
if "$meta_executable" \
    --input_a "/nonexistent/file_a.bed" \
    --input_b "$meta_temp_dir/input_b.bed" \
    --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "✗ Should have failed with non-existent input_a file"
    exit 1
else
    log "✓ Properly handles non-existent input_a files"
fi

if "$meta_executable" \
    --input_a "$meta_temp_dir/input_a.bed" \
    --input_b "/nonexistent/file_b.bed" \
    --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "✗ Should have failed with non-existent input_b file"
    exit 1
else
    log "✓ Properly handles non-existent input_b files"
fi

####################################################################################################

log "TEST 10: Empty input handling"

# Create empty input files
touch "$meta_temp_dir/empty_a.bed"
touch "$meta_temp_dir/empty_b.bed"

# Empty A file should produce empty output
"$meta_executable" \
    --input_a "$meta_temp_dir/empty_a.bed" \
    --input_b "$meta_temp_dir/input_b.bed" \
    --output "$meta_temp_dir/empty_a_output.bed"

check_file_exists "$meta_temp_dir/empty_a_output.bed" "empty A input test output"

if [ ! -s "$meta_temp_dir/empty_a_output.bed" ]; then
    log "✓ Empty input A produces empty output"
else
    log "✗ Empty input A handling failed"
    cat "$meta_temp_dir/empty_a_output.bed"
    exit 1
fi

# Empty B file should preserve all A intervals
"$meta_executable" \
    --input_a "$meta_temp_dir/input_a.bed" \
    --input_b "$meta_temp_dir/empty_b.bed" \
    --output "$meta_temp_dir/empty_b_output.bed"

check_file_exists "$meta_temp_dir/empty_b_output.bed" "empty B input test output"

# Output should match input A exactly
input_lines=$(wc -l < "$meta_temp_dir/input_a.bed")
output_lines=$(wc -l < "$meta_temp_dir/empty_b_output.bed")

if [ "$input_lines" -eq "$output_lines" ]; then
    log "✓ Empty input B preserves all A intervals ($output_lines lines)"
else
    log "✗ Empty input B handling failed: expected $input_lines, got $output_lines"
    exit 1
fi

####################################################################################################

log "TEST 11: No overlap scenario"

# Create test data with no overlapping intervals
cat > "$meta_temp_dir/no_overlap_a.bed" << 'EOF'
chr1	100	200	feature1	100	+
chr1	300	400	feature2	200	+
chr1	500	600	feature3	300	+
EOF

cat > "$meta_temp_dir/no_overlap_b.bed" << 'EOF'
chr1	50	90	distant1	50	.
chr1	210	290	distant2	75	.
chr1	410	490	distant3	100	.
chr1	650	750	distant4	150	.
EOF

"$meta_executable" \
    --input_a "$meta_temp_dir/no_overlap_a.bed" \
    --input_b "$meta_temp_dir/no_overlap_b.bed" \
    --output "$meta_temp_dir/no_overlap_output.bed"

check_file_exists "$meta_temp_dir/no_overlap_output.bed" "no overlap output"

# All intervals should be preserved since there are no overlaps
input_lines=$(wc -l < "$meta_temp_dir/no_overlap_a.bed")
output_lines=$(wc -l < "$meta_temp_dir/no_overlap_output.bed")

if [ "$input_lines" -eq "$output_lines" ]; then
    log "✓ No overlap scenario preserves all intervals ($output_lines lines)"
else
    log "✗ No overlap handling failed: expected $input_lines, got $output_lines"
    exit 1
fi

####################################################################################################

log "TEST 12: Multiple overlapping intervals from B"

# Create scenario where one A interval overlaps with multiple B intervals
cat > "$meta_temp_dir/multi_overlap_a.bed" << 'EOF'
chr1	100	500	big_feature	1000	+
EOF

cat > "$meta_temp_dir/multi_overlap_b.bed" << 'EOF'
chr1	150	200	repeat1	50	.
chr1	250	300	repeat2	75	.
chr1	350	400	repeat3	100	.
EOF

"$meta_executable" \
    --input_a "$meta_temp_dir/multi_overlap_a.bed" \
    --input_b "$meta_temp_dir/multi_overlap_b.bed" \
    --output "$meta_temp_dir/multi_overlap_output.bed"

check_file_exists "$meta_temp_dir/multi_overlap_output.bed" "multiple overlap output"
check_file_not_empty "$meta_temp_dir/multi_overlap_output.bed" "multiple overlap output"

# Should create multiple intervals: 100-150, 200-250, 300-350, 400-500
expected_intervals=4
actual_intervals=$(wc -l < "$meta_temp_dir/multi_overlap_output.bed")

if [ "$actual_intervals" -eq "$expected_intervals" ]; then
    log "✓ Multiple overlapping B intervals create correct number of output intervals ($actual_intervals)"
else
    log "✗ Multiple overlap handling failed: expected $expected_intervals, got $actual_intervals"
    cat "$meta_temp_dir/multi_overlap_output.bed"
    exit 1
fi

####################################################################################################

log "All tests completed successfully!"
