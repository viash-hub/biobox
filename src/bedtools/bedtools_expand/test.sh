#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_expand"

# Create test data
log "Creating test data..."

# Create simple test file with comma-separated values in one column
cat > "$meta_temp_dir/simple.bed" << 'EOF'
chr1	100	200	1,2,3
chr1	300	400	4,5,6
chr2	500	600	7,8
EOF

# Create test file with comma-separated values in multiple columns
cat > "$meta_temp_dir/multi_column.bed" << 'EOF'
chr1	10	20	1,2,3	10,20,30
chr1	40	50	4,5,6	40,50,60
chr2	70	80	7,8,9	70,80,90
EOF

# Create BED file with single values (no expansion needed)
cat > "$meta_temp_dir/no_expansion.bed" << 'EOF'
chr1	100	200	single_value
chr2	300	400	another_value
EOF

# Create file with unequal comma-separated lists (should be handled gracefully)
cat > "$meta_temp_dir/unequal.bed" << 'EOF'
chr1	100	200	1,2,3	10,20
chr1	300	400	4,5	40,50,60
EOF

# Test 1: Basic single column expansion
log "Starting TEST 1: Basic single column expansion"
"$meta_executable" \
    --input "$meta_temp_dir/simple.bed" \
    --columns "4" \
    --output "$meta_temp_dir/output1.bed"

check_file_exists "$meta_temp_dir/output1.bed" "single column expansion output"
check_file_not_empty "$meta_temp_dir/output1.bed" "single column expansion output"
check_file_line_count "$meta_temp_dir/output1.bed" 8 "single column expansion line count"

# Check that expansion worked correctly
check_file_contains "$meta_temp_dir/output1.bed" "chr1	100	200	1" "first expanded value"
check_file_contains "$meta_temp_dir/output1.bed" "chr1	100	200	2" "second expanded value"
check_file_contains "$meta_temp_dir/output1.bed" "chr1	100	200	3" "third expanded value"
check_file_contains "$meta_temp_dir/output1.bed" "chr2	500	600	7" "chr2 first value"
check_file_contains "$meta_temp_dir/output1.bed" "chr2	500	600	8" "chr2 second value"
log "✅ TEST 1 completed successfully"

# Test 2: Multi-column expansion
log "Starting TEST 2: Multi-column expansion"
"$meta_executable" \
    --input "$meta_temp_dir/multi_column.bed" \
    --columns "4,5" \
    --output "$meta_temp_dir/output2.bed"

check_file_exists "$meta_temp_dir/output2.bed" "multi-column expansion output"
check_file_not_empty "$meta_temp_dir/output2.bed" "multi-column expansion output"
check_file_line_count "$meta_temp_dir/output2.bed" 9 "multi-column expansion line count"

# Check that paired expansion worked correctly
check_file_contains "$meta_temp_dir/output2.bed" "chr1	10	20	1	10" "first paired expansion"
check_file_contains "$meta_temp_dir/output2.bed" "chr1	10	20	2	20" "second paired expansion"
check_file_contains "$meta_temp_dir/output2.bed" "chr1	10	20	3	30" "third paired expansion"
log "✅ TEST 2 completed successfully"

# Test 3: No expansion needed (single values)
log "Starting TEST 3: Single values (no expansion needed)"
"$meta_executable" \
    --input "$meta_temp_dir/no_expansion.bed" \
    --columns "4" \
    --output "$meta_temp_dir/output3.bed"

check_file_exists "$meta_temp_dir/output3.bed" "no expansion output"
check_file_not_empty "$meta_temp_dir/output3.bed" "no expansion output"
check_file_line_count "$meta_temp_dir/output3.bed" 2 "no expansion line count"

# Should be identical to input since no comma-separated values
check_file_contains "$meta_temp_dir/output3.bed" "single_value" "single value preserved"
check_file_contains "$meta_temp_dir/output3.bed" "another_value" "another value preserved"
log "✅ TEST 3 completed successfully"

# Test 4: Different column positions
log "Starting TEST 4: Different column positions"
"$meta_executable" \
    --input "$meta_temp_dir/multi_column.bed" \
    --columns "5" \
    --output "$meta_temp_dir/output4.bed"

check_file_exists "$meta_temp_dir/output4.bed" "column 5 expansion output"
check_file_not_empty "$meta_temp_dir/output4.bed" "column 5 expansion output"
check_file_line_count "$meta_temp_dir/output4.bed" 9 "column 5 expansion line count"

# Check that only column 5 was expanded, column 4 remains comma-separated
check_file_contains "$meta_temp_dir/output4.bed" "chr1	10	20	1,2,3	10" "column 4 not expanded"
check_file_contains "$meta_temp_dir/output4.bed" "chr1	10	20	1,2,3	20" "column 5 expanded"
log "✅ TEST 4 completed successfully"

# Test 5: Large expansion test
log "Starting TEST 5: Large expansion test"
# Create file with more comma-separated values
cat > "$meta_temp_dir/large.bed" << 'EOF'
chr1	100	200	1,2,3,4,5,6,7,8,9,10
EOF

"$meta_executable" \
    --input "$meta_temp_dir/large.bed" \
    --columns "4" \
    --output "$meta_temp_dir/output5.bed"

check_file_exists "$meta_temp_dir/output5.bed" "large expansion output"
check_file_not_empty "$meta_temp_dir/output5.bed" "large expansion output"
check_file_line_count "$meta_temp_dir/output5.bed" 10 "large expansion line count"

# Check that all values are expanded
for i in {1..10}; do
    if ! grep -q "chr1	100	200	$i$" "$meta_temp_dir/output5.bed"; then
        log_error "Expected value $i not found in large expansion"
        exit 1
    fi
done
log "✅ TEST 5 completed successfully"

log "🎉 All bedtools_expand tests completed successfully!"
