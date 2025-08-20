#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_sample"

####################################################################################################

log "Creating test data..."

# Create test BED file with multiple records for sampling
create_test_bed "$meta_temp_dir/input.bed" 1000
check_file_exists "$meta_temp_dir/input.bed" "input BED file"

# Verify test data was created correctly
input_lines=$(wc -l < "$meta_temp_dir/input.bed")
log "Created test BED file with $input_lines records"

####################################################################################################

log "TEST 1: Basic sampling functionality"

# Test basic sampling with specified number (smaller than input)
"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --output "$meta_temp_dir/basic_sample.bed" \
    --number 100

check_file_exists "$meta_temp_dir/basic_sample.bed" "basic sample output"

# Check that output has records (should be min of input size or requested amount)
output_lines=$(wc -l < "$meta_temp_dir/basic_sample.bed")
log "✓ Basic sampling produced $output_lines records"

if [ "$output_lines" -eq 100 ]; then
    log "✓ Output record count is as expected (100)"
else
    log "✗ Expected 100 records, got $output_lines"
    exit 1
fi

####################################################################################################

log "TEST 2: Specific number sampling"

# Test sampling with specific number
sample_count=50
"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --output "$meta_temp_dir/numbered_sample.bed" \
    --number "$sample_count"

check_file_exists "$meta_temp_dir/numbered_sample.bed" "numbered sample output"

# Check that output has exactly the requested number of records
output_lines=$(wc -l < "$meta_temp_dir/numbered_sample.bed")
log "Requested $sample_count records, got $output_lines records"

if [ "$output_lines" -eq "$sample_count" ]; then
    log "✓ Exact number sampling works correctly"
else
    log "✗ Did not get expected number of records"
    exit 1
fi

####################################################################################################

log "TEST 3: Reproducible sampling with seed"

# Test that same seed produces same results
seed_value=123
"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --output "$meta_temp_dir/seed_sample1.bed" \
    --number 25 \
    --seed "$seed_value"

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --output "$meta_temp_dir/seed_sample2.bed" \
    --number 25 \
    --seed "$seed_value"

check_file_exists "$meta_temp_dir/seed_sample1.bed" "first seeded sample"
check_file_exists "$meta_temp_dir/seed_sample2.bed" "second seeded sample"

# Check that both outputs are identical
if diff -q "$meta_temp_dir/seed_sample1.bed" "$meta_temp_dir/seed_sample2.bed" >/dev/null; then
    log "✓ Seeded sampling produces reproducible results"
else
    log "✗ Seeded sampling is not reproducible"
    exit 1
fi

####################################################################################################

log "TEST 4: Over-sampling (request more records than available)"

# Test when requested sample size exceeds input size
small_input="$meta_temp_dir/small_input.bed"
create_test_bed "$small_input" 10

# bedtools sample returns an error when requesting more records than available
if "$meta_executable" \
    --input "$small_input" \
    --output "$meta_temp_dir/oversample.bed" \
    --number 100 2>/dev/null; then
    log "✗ Should have failed when requesting more records than available"
    exit 1
else
    log "✓ Correctly handles over-sampling by returning error (expected behavior)"
fi

####################################################################################################

log "TEST 5: Different random seeds produce different results"

# Test that different seeds produce different samples
"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --output "$meta_temp_dir/diff_seed1.bed" \
    --number 30 \
    --seed 456

"$meta_executable" \
    --input "$meta_temp_dir/input.bed" \
    --output "$meta_temp_dir/diff_seed2.bed" \
    --number 30 \
    --seed 789

check_file_exists "$meta_temp_dir/diff_seed1.bed" "first different seed sample"
check_file_exists "$meta_temp_dir/diff_seed2.bed" "second different seed sample"

# Different seeds should produce different samples (with high probability)
if ! diff -q "$meta_temp_dir/diff_seed1.bed" "$meta_temp_dir/diff_seed2.bed" >/dev/null; then
    log "✓ Different seeds produce different samples"
else
    log "⚠ Different seeds produced identical samples (possible but unlikely)"
fi

####################################################################################################

log "TEST 6: Header inclusion option"

# Create test data with header-like content
cat > "$meta_temp_dir/with_header.bed" << 'EOF'
# This is a header line
# Another header line
chr1	100	200	feature1	100	+
chr1	300	400	feature2	200	-
chr2	150	250	feature3	150	+
EOF

"$meta_executable" \
    --input "$meta_temp_dir/with_header.bed" \
    --output "$meta_temp_dir/header_test.bed" \
    --number 2 \
    --include_header \
    --seed 999

check_file_exists "$meta_temp_dir/header_test.bed" "header inclusion test"

# Check if output preserves the format (may include comments)
output_lines=$(wc -l < "$meta_temp_dir/header_test.bed")
log "✓ Header inclusion test completed with $output_lines lines"

####################################################################################################

log "TEST 7: Parameter validation"

# Test required parameter validation
log "Testing parameter validation"

if "$meta_executable" --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "✗ Should have failed without --input parameter"
    exit 1
else
    log "✓ Correctly requires --input parameter"
fi

if "$meta_executable" --input "$meta_temp_dir/input.bed" 2>/dev/null; then
    log "✗ Should have failed without --output parameter"  
    exit 1
else
    log "✓ Correctly requires --output parameter"
fi

####################################################################################################

log "TEST 8: File format validation"

# Test with non-existent input file
if "$meta_executable" \
    --input "/nonexistent/file.bed" \
    --output "$meta_temp_dir/test.bed" 2>/dev/null; then
    log "Should have failed with non-existent input file"
else
    log "✓ Properly handles non-existent input files"
fi

####################################################################################################

log "All tests completed successfully!"
