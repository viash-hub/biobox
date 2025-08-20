#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_split"

####################################################################################################

log "Creating test data..."

# Create test BED file with multiple intervals for splitting
cat > "$meta_temp_dir/input.bed" << 'EOF'
chr1	0	100	feature1	100	+
chr1	200	300	feature2	200	+
chr1	400	500	feature3	300	+
chr1	600	700	feature4	400	+
chr1	800	900	feature5	500	+
chr1	1000	1100	feature6	600	+
chr2	0	50	feature7	700	-
chr2	100	150	feature8	800	-
chr2	200	250	feature9	900	-
chr2	300	350	feature10	1000	-
EOF

check_file_exists "$meta_temp_dir/input.bed" "input BED file"

####################################################################################################

log "TEST 1: Basic split into 2 files"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --number 2 \
    --output_dir "$meta_temp_dir" \
    --prefix "split_2"

# bedtools split creates files with format: prefix.00001.bed, prefix.00002.bed, etc.
check_file_exists "$meta_temp_dir/split_2.00001.bed" "first split file"
check_file_exists "$meta_temp_dir/split_2.00002.bed" "second split file"
check_file_not_empty "$meta_temp_dir/split_2.00001.bed" "first split file"
check_file_not_empty "$meta_temp_dir/split_2.00002.bed" "second split file"

# Count total lines to ensure all records are preserved
total_input=$(wc -l < "$meta_temp_dir/input.bed")
total_output=$(($(wc -l < "$meta_temp_dir/split_2.00001.bed") + $(wc -l < "$meta_temp_dir/split_2.00002.bed")))

if [ "$total_input" -eq "$total_output" ]; then
    log "✓ All records preserved in split files ($total_output)"
else
    log "✗ Record count mismatch: input=$total_input, output=$total_output"
    exit 1
fi

####################################################################################################

log "TEST 2: Split into 3 files with size algorithm"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --number 3 \
    --output_dir "$meta_temp_dir" \
    --prefix "split_3_size" \
    --algorithm "size"

check_file_exists "$meta_temp_dir/split_3_size.00001.bed" "first file (size algorithm)"
check_file_exists "$meta_temp_dir/split_3_size.00002.bed" "second file (size algorithm)"
check_file_exists "$meta_temp_dir/split_3_size.00003.bed" "third file (size algorithm)"

# Verify all files have content
for i in {1..3}; do
    file_num=$(printf "%05d" $i)
    check_file_not_empty "$meta_temp_dir/split_3_size.$file_num.bed" "file $i (size algorithm)"
done

# Count total lines
total_output=$(($(wc -l < "$meta_temp_dir/split_3_size.00001.bed") + \
                $(wc -l < "$meta_temp_dir/split_3_size.00002.bed") + \
                $(wc -l < "$meta_temp_dir/split_3_size.00003.bed")))

if [ "$total_input" -eq "$total_output" ]; then
    log "✓ All records preserved with size algorithm ($total_output)"
else
    log "✗ Record count mismatch with size algorithm: input=$total_input, output=$total_output"
    exit 1
fi

####################################################################################################

log "TEST 3: Split into 3 files with simple algorithm"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --number 3 \
    --output_dir "$meta_temp_dir" \
    --prefix "split_3_simple" \
    --algorithm "simple"

check_file_exists "$meta_temp_dir/split_3_simple.00001.bed" "first file (simple algorithm)"
check_file_exists "$meta_temp_dir/split_3_simple.00002.bed" "second file (simple algorithm)"
check_file_exists "$meta_temp_dir/split_3_simple.00003.bed" "third file (simple algorithm)"

# With simple algorithm, files should have approximately equal number of records
file1_lines=$(wc -l < "$meta_temp_dir/split_3_simple.00001.bed")
file2_lines=$(wc -l < "$meta_temp_dir/split_3_simple.00002.bed")
file3_lines=$(wc -l < "$meta_temp_dir/split_3_simple.00003.bed")

total_simple=$((file1_lines + file2_lines + file3_lines))

if [ "$total_input" -eq "$total_simple" ]; then
    log "✓ All records preserved with simple algorithm ($total_simple)"
    log "  File 1: $file1_lines lines, File 2: $file2_lines lines, File 3: $file3_lines lines"
else
    log "✗ Record count mismatch with simple algorithm: input=$total_input, output=$total_simple"
    exit 1
fi

# Check that files have roughly equal numbers of records (within 1-2 of each other)
expected_per_file=$((total_input / 3))
for lines in $file1_lines $file2_lines $file3_lines; do
    diff=$((lines - expected_per_file))
    if [ ${diff#-} -le 2 ]; then  # abs(diff) <= 2
        continue
    else
        log "✗ Simple algorithm distribution not balanced: expected ~$expected_per_file, got $lines"
        exit 1
    fi
done
log "✓ Simple algorithm produces balanced distribution"

####################################################################################################

log "TEST 4: Split with default prefix"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --number 2 \
    --output_dir "$meta_temp_dir"

# When no prefix specified, bedtools uses "_split" as default prefix
check_file_exists "$meta_temp_dir/_split.00001.bed" "first file (default prefix)"
check_file_exists "$meta_temp_dir/_split.00002.bed" "second file (default prefix)"

####################################################################################################

log "TEST 5: Parameter validation"

# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" \
    --number 2 \
    --output_dir "$meta_temp_dir" 2>/dev/null; then
    log "✗ Should have failed without --input parameter"
    exit 1
else
    log "✓ Correctly requires --input parameter"
fi

if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --output_dir "$meta_temp_dir" 2>/dev/null; then
    log "✗ Should have failed without --number parameter"
    exit 1
else
    log "✓ Correctly requires --number parameter"
fi

if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --number 2 2>/dev/null; then
    log "✗ Should have failed without --output_dir parameter"
    exit 1
else
    log "✓ Correctly requires --output_dir parameter"
fi

####################################################################################################

log "TEST 6: Algorithm parameter validation"

# Test invalid algorithm
if "$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --number 2 \
    --output_dir "$meta_temp_dir" \
    --prefix "invalid_test" \
    --algorithm "invalid" 2>/dev/null; then
    log "✗ Should have failed with invalid algorithm"
    exit 1
else
    log "✓ Correctly rejects invalid algorithm"
fi

####################################################################################################

log "TEST 7: File validation"

# Test with non-existent files
if "$meta_executable" \
    --input "/nonexistent/file.bed" \
    --number 2 \
    --output_dir "$meta_temp_dir" \
    --prefix "test" 2>/dev/null; then
    log "✗ Should have failed with non-existent input file"
    exit 1
else
    log "✓ Properly handles non-existent input files"
fi

####################################################################################################

log "TEST 8: Empty input handling"

# Create empty input file
touch "$meta_temp_dir/empty.bed"

# bedtools split will give a warning but should not fail
"$meta_executable" \
    --input "$meta_temp_dir/empty.bed" \
    --number 2 \
    --output_dir "$meta_temp_dir" \
    --prefix "empty_test"

# For empty input, bedtools split doesn't create any output files
# This is expected behavior - we just verify the command doesn't crash
if [ ! -f "$meta_temp_dir/empty_test.00001.bed" ]; then
    log "✓ Empty input correctly produces no output files (expected behavior)"
else
    # If files were created, they should be empty
    if [ ! -s "$meta_temp_dir/empty_test.00001.bed" ] && [ ! -s "$meta_temp_dir/empty_test.00002.bed" ]; then
        log "✓ Empty input produces empty split files"
    else
        log "✗ Empty input handling failed"
        exit 1
    fi
fi

####################################################################################################

log "TEST 9: Single record split"

# Create file with single record
cat > "$meta_temp_dir/single.bed" << 'EOF'
chr1	100	200	single_feature	100	+
EOF

"$meta_executable" \
    --input "$meta_temp_dir/single.bed" \
    --number 3 \
    --output_dir "$meta_temp_dir" \
    --prefix "single_test"

# When input has fewer records than requested files, bedtools split only creates
# as many files as needed. With 1 record and 3 files requested, only 1 file is created.
check_file_exists "$meta_temp_dir/single_test.00001.bed" "single record split file 1"
check_file_not_empty "$meta_temp_dir/single_test.00001.bed" "single record split file 1"

# Verify the single record is in the first file
if [ "$(wc -l < "$meta_temp_dir/single_test.00001.bed")" -eq 1 ]; then
    log "✓ Single record correctly placed in first split file"
else
    log "✗ Single record split failed"
    exit 1
fi

# Check that no additional files were created (this is expected behavior)
if [ ! -f "$meta_temp_dir/single_test.00002.bed" ] && [ ! -f "$meta_temp_dir/single_test.00003.bed" ]; then
    log "✓ No unnecessary empty files created for single record"
else
    log "✓ Additional files created (may be empty, which is also acceptable)"
fi

####################################################################################################

log "TEST 10: Large number split test"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --number 5 \
    --output_dir "$meta_temp_dir" \
    --prefix "split_5"

# Should create 5 files
for i in {1..5}; do
    file_num=$(printf "%05d" $i)
    check_file_exists "$meta_temp_dir/split_5.$file_num.bed" "split file $i"
done

# Count total records
total_split5=0
for i in {1..5}; do
    file_num=$(printf "%05d" $i)
    if [ -s "$meta_temp_dir/split_5.$file_num.bed" ]; then
        lines=$(wc -l < "$meta_temp_dir/split_5.$file_num.bed")
        total_split5=$((total_split5 + lines))
    fi
done

if [ "$total_input" -eq "$total_split5" ]; then
    log "✓ All records preserved when splitting into 5 files ($total_split5)"
else
    log "✗ Record count mismatch for 5-way split: input=$total_input, output=$total_split5"
    exit 1
fi

####################################################################################################

log "All tests completed successfully!"
