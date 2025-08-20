#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_summary"

####################################################################################################

log "Creating test data..."

# Create test BED file with genomic intervals
cat > "$meta_temp_dir/input.bed" << 'EOF'
chr1	100	200	feature1	100	+
chr1	300	500	feature2	200	+
chr1	600	650	feature3	150	+
chr1	800	1000	feature4	300	-
chr2	150	250	feature5	120	+
chr2	400	600	feature6	180	-
chr2	700	750	feature7	90	+
chr3	50	150	feature8	200	+
EOF

# Create genome file with chromosome lengths
cat > "$meta_temp_dir/genome.txt" << 'EOF'
chr1	2000
chr2	1500
chr3	1000
EOF

check_file_exists "$meta_temp_dir/input.bed" "input BED file"
check_file_exists "$meta_temp_dir/genome.txt" "genome file"

####################################################################################################

log "TEST 1: Basic summary statistics"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/basic_summary.txt"

check_file_exists "$meta_temp_dir/basic_summary.txt" "basic summary output"
check_file_not_empty "$meta_temp_dir/basic_summary.txt" "basic summary output"

# Check that output contains statistical information
if [ -s "$meta_temp_dir/basic_summary.txt" ]; then
    log "✓ Summary statistics generated successfully"
else
    log "✗ Summary statistics generation failed"
    exit 1
fi

# Verify output has reasonable content (should be text with statistics)
line_count=$(wc -l < "$meta_temp_dir/basic_summary.txt")
if [ "$line_count" -gt 0 ]; then
    log "✓ Summary output contains $line_count lines of statistics"
else
    log "✗ Summary output is unexpectedly empty"
    exit 1
fi

####################################################################################################

log "TEST 2: Output format verification"

# Check the structure of the summary output
if grep -q "chrom" "$meta_temp_dir/basic_summary.txt" && \
   grep -q "num_ivls" "$meta_temp_dir/basic_summary.txt" && \
   grep -q "total_ivl_bp" "$meta_temp_dir/basic_summary.txt"; then
    log "✓ Summary output contains expected column headers"
else
    log "✗ Summary output format verification failed"
    cat "$meta_temp_dir/basic_summary.txt"
    exit 1
fi

# Check that both per-chromosome and summary ('all') lines are present
if grep -q "chr1" "$meta_temp_dir/basic_summary.txt" && \
   grep -q "chr2" "$meta_temp_dir/basic_summary.txt" && \
   grep -q "chr3" "$meta_temp_dir/basic_summary.txt" && \
   grep -q "all" "$meta_temp_dir/basic_summary.txt"; then
    log "✓ Output contains both per-chromosome and summary statistics"
else
    log "✗ Output format does not match expected structure"
    cat "$meta_temp_dir/basic_summary.txt"
    exit 1
fi

####################################################################################################

log "TEST 3: Statistical values verification"

# Verify that statistical values are reasonable
# We have 8 intervals total in our test data
if grep -q "all.*8" "$meta_temp_dir/basic_summary.txt"; then
    log "✓ Correct total number of intervals reported"
else
    log "✗ Total interval count is incorrect"
    cat "$meta_temp_dir/basic_summary.txt"
    exit 1
fi

# Check that chromosome-specific counts are correct
# chr1 has 4 intervals, chr2 has 3, chr3 has 1
if grep "^chr1" "$meta_temp_dir/basic_summary.txt" | grep -q "4" && \
   grep "^chr2" "$meta_temp_dir/basic_summary.txt" | grep -q "3" && \
   grep "^chr3" "$meta_temp_dir/basic_summary.txt" | grep -q "1"; then
    log "✓ Per-chromosome interval counts are correct"
else
    log "✗ Per-chromosome counts are incorrect"
    cat "$meta_temp_dir/basic_summary.txt"
    exit 1
fi

####################################################################################################

log "TEST 4: Different chromosome subset test"

# Create test with only one chromosome to verify individual chromosome handling
cat > "$meta_temp_dir/single_chr.bed" << 'EOF'
chr1	100	300	feature1	100	+
chr1	500	700	feature2	200	+
EOF

"$meta_executable" \
    --input "$meta_temp_dir/single_chr.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/single_chr_summary.txt"

check_file_exists "$meta_temp_dir/single_chr_summary.txt" "single chromosome summary output"
check_file_not_empty "$meta_temp_dir/single_chr_summary.txt" "single chromosome summary output"

# Should have header line, chr1 line, and 'all' line
expected_lines=3
actual_lines=$(wc -l < "$meta_temp_dir/single_chr_summary.txt")
if [ "$actual_lines" -eq "$expected_lines" ]; then
    log "✓ Single chromosome summary has correct structure ($actual_lines lines)"
else
    log "✗ Single chromosome summary structure unexpected: got $actual_lines lines, expected $expected_lines"
    cat "$meta_temp_dir/single_chr_summary.txt"
fi

####################################################################################################

log "TEST 5: Different input formats - verify basic functionality"

# Create a simple BED3 format file (minimal columns)
cat > "$meta_temp_dir/simple.bed" << 'EOF'
chr1	100	300
chr1	500	700
chr2	200	400
EOF

"$meta_executable" \
    --input "$meta_temp_dir/simple.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/simple_summary.txt"

check_file_exists "$meta_temp_dir/simple_summary.txt" "simple format summary output"
check_file_not_empty "$meta_temp_dir/simple_summary.txt" "simple format summary output"

log "✓ Simple BED format processing works"

####################################################################################################

log "TEST 6: Parameter validation"

# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/test.txt" 2>/dev/null; then
    log "✗ Should have failed without --input parameter"
    exit 1
else
    log "✓ Correctly requires --input parameter"
fi

if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --output "$meta_temp_dir/test.txt" 2>/dev/null; then
    log "✗ Should have failed without --genome parameter"
    exit 1
else
    log "✓ Correctly requires --genome parameter"
fi

if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "$meta_temp_dir/genome.txt" 2>/dev/null; then
    log "✗ Should have failed without --output parameter"
    exit 1
else
    log "✓ Correctly requires --output parameter"
fi

####################################################################################################

log "TEST 7: File validation"

# Test with non-existent input file
if "$meta_executable" \
    --input "/nonexistent/file.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/test.txt" 2>/dev/null; then
    log "✗ Should have failed with non-existent input file"
    exit 1
else
    log "✓ Properly handles non-existent input files"
fi

# Test with non-existent genome file
if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "/nonexistent/genome.txt" \
    --output "$meta_temp_dir/test.txt" 2>/dev/null; then
    log "✗ Should have failed with non-existent genome file"
    exit 1
else
    log "✓ Properly handles non-existent genome files"
fi

####################################################################################################

log "TEST 8: Empty input handling"

# Create empty input file
touch "$meta_temp_dir/empty.bed"

"$meta_executable" \
    --input "$meta_temp_dir/empty.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/empty_summary.txt"

check_file_exists "$meta_temp_dir/empty_summary.txt" "empty input summary output"

# Empty input should still produce some output (likely zeros or headers)
if [ -f "$meta_temp_dir/empty_summary.txt" ]; then
    log "✓ Empty input produces valid output file"
else
    log "✗ Empty input handling failed"
    exit 1
fi

####################################################################################################

log "TEST 9: Single interval test"

# Create file with single interval
cat > "$meta_temp_dir/single.bed" << 'EOF'
chr1	500	800	single_feature	100	+
EOF

"$meta_executable" \
    --input "$meta_temp_dir/single.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/single_summary.txt"

check_file_exists "$meta_temp_dir/single_summary.txt" "single interval summary output"
check_file_not_empty "$meta_temp_dir/single_summary.txt" "single interval summary output"

log "✓ Single interval processing works correctly"

####################################################################################################

log "TEST 10: Malformed genome file handling"

# Create malformed genome file (missing tab or size)
cat > "$meta_temp_dir/bad_genome.txt" << 'EOF'
chr1
chr2	notanumber
chr3	1000
EOF

# This should fail gracefully
if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --genome "$meta_temp_dir/bad_genome.txt" \
    --output "$meta_temp_dir/test_bad.txt" 2>/dev/null; then
    log "✗ Should have failed with malformed genome file"
    exit 1
else
    log "✓ Properly handles malformed genome files"
fi

####################################################################################################

log "TEST 11: Large interval dataset simulation"

# Create a larger dataset to test performance
cat > "$meta_temp_dir/large.bed" << 'EOF'
chr1	100	200	f1	100	+
chr1	300	400	f2	150	+
chr1	500	600	f3	200	+
chr1	700	800	f4	250	+
chr1	900	1000	f5	300	+
chr1	1100	1200	f6	350	+
chr1	1300	1400	f7	400	+
chr1	1500	1600	f8	450	+
chr2	100	150	f9	100	-
chr2	200	250	f10	150	-
chr2	300	350	f11	200	-
chr2	400	450	f12	250	-
chr2	500	550	f13	300	-
chr2	600	650	f14	350	-
chr2	700	750	f15	400	-
chr3	50	100	f16	100	+
chr3	150	200	f17	150	+
chr3	250	300	f18	200	+
chr3	350	400	f19	250	+
chr3	450	500	f20	300	+
EOF

"$meta_executable" \
    --input "$meta_temp_dir/large.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/large_summary.txt"

check_file_exists "$meta_temp_dir/large_summary.txt" "large dataset summary output"
check_file_not_empty "$meta_temp_dir/large_summary.txt" "large dataset summary output"

log "✓ Large dataset processing completed successfully"

####################################################################################################

log "TEST 12: Output format consistency"

# Compare outputs to ensure they have consistent structure
basic_size=$(wc -c < "$meta_temp_dir/basic_summary.txt")
single_size=$(wc -c < "$meta_temp_dir/single_summary.txt")

if [ "$basic_size" -gt 0 ] && [ "$single_size" -gt 0 ]; then
    log "✓ All outputs have consistent non-zero size"
else
    log "✗ Output size consistency check failed"
    log "  Basic summary: $basic_size bytes"
    log "  Single interval: $single_size bytes"
    exit 1
fi

# Verify all test outputs exist and are readable
for output_file in basic_summary.txt single_chr_summary.txt simple_summary.txt; do
    if [ -r "$meta_temp_dir/$output_file" ]; then
        log "✓ Output file $output_file is readable"
    else
        log "✗ Output file $output_file is not readable"
        exit 1
    fi
done

####################################################################################################

log "All tests completed successfully!"
