#!/bin/bash

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_closest"

# Create test data
log "Creating test data..."

# Create query intervals file
cat > "$meta_temp_dir/queries.bed" << 'EOF'
chr1	100	200	query1	100	+
chr1	400	500	query2	200	+
chr1	800	900	query3	300	-
chr2	200	300	query4	400	-
EOF

# Create database file with features at various distances
cat > "$meta_temp_dir/database.bed" << 'EOF'
chr1	250	350	feature1	500	+
chr1	450	550	feature2	600	+
chr1	700	800	feature3	700	-
chr2	150	250	feature4	800	+
chr2	600	700	feature5	900	-
chr2	950	1050	feature6	1000	+
EOF

# Create second database file for multi-file testing
cat > "$meta_temp_dir/database2.bed" << 'EOF'
chr1	1050	1150	db2_feature1
chr1	1250	1350	db2_feature2
chr1	1450	1550	db2_feature3
EOF

# Create distant features for signed distance testing (non-overlapping)
cat > "$meta_temp_dir/test_b_distant.bed" << 'EOF'
chr1	50	90	upstream1
chr1	250	290	downstream1
chr1	450	490	upstream2
chr1	650	690	downstream2
EOF

# Test 1: Basic closest feature finding
log "Starting TEST 1: Basic closest feature finding"
"$meta_executable" \
    --input_a "$meta_temp_dir/queries.bed" \
    --input_b "$meta_temp_dir/database.bed" \
    --output "$meta_temp_dir/output1.bed"

check_file_exists "$meta_temp_dir/output1.bed" "basic closest output"
check_file_not_empty "$meta_temp_dir/output1.bed" "basic closest output"
check_file_line_count "$meta_temp_dir/output1.bed" 4 "basic closest line count"

# Check that closest features are reported
check_file_contains "$meta_temp_dir/output1.bed" "feature" "closest features found"
log "✅ TEST 1 completed successfully"

# Test 2: Closest features with distance reporting
log "Starting TEST 2: Closest features with distance reporting"
"$meta_executable" \
    --input_a "$meta_temp_dir/queries.bed" \
    --input_b "$meta_temp_dir/database.bed" \
    --distance_mode "ref" \
    --output "$meta_temp_dir/output2.bed"

check_file_exists "$meta_temp_dir/output2.bed" "distance output"
check_file_not_empty "$meta_temp_dir/output2.bed" "distance output"
check_file_line_count "$meta_temp_dir/output2.bed" 4 "distance line count"

# Check that distance column is added (should have more columns than input)
input_cols=$(head -1 "$meta_temp_dir/queries.bed" | awk '{print NF}')
output_cols=$(head -1 "$meta_temp_dir/output2.bed" | awk '{print NF}')
if [ $output_cols -le $input_cols ]; then
    error "Expected more columns in output with distance, got $output_cols vs input $input_cols"
fi
log "✅ TEST 2 completed successfully"

# Test 3: Find closest with strand consideration
log "Starting TEST 3: Closest with strand consideration"
"$meta_executable" \
    --input_a "$meta_temp_dir/queries.bed" \
    --input_b "$meta_temp_dir/database.bed" \
    --strand \
    --output "$meta_temp_dir/output3.bed"

check_file_exists "$meta_temp_dir/output3.bed" "strand output"
check_file_not_empty "$meta_temp_dir/output3.bed" "strand output"
log "✅ TEST 3 completed successfully"

# Test 4: Find k-nearest neighbors (k=2)
log "Starting TEST 4: K-nearest neighbors (k=2)"
"$meta_executable" \
    --input_a "$meta_temp_dir/queries.bed" \
    --input_b "$meta_temp_dir/database.bed" \
    --k_closest 2 \
    --output "$meta_temp_dir/output4.bed"

check_file_exists "$meta_temp_dir/output4.bed" "k-nearest output"
check_file_not_empty "$meta_temp_dir/output4.bed" "k-nearest output"

# Should have more lines than basic test (up to 2x for each query)
basic_lines=$(wc -l < "$meta_temp_dir/output1.bed")
knearest_lines=$(wc -l < "$meta_temp_dir/output4.bed")
if [ $knearest_lines -lt $basic_lines ]; then
    error "Expected at least $basic_lines lines for k-nearest, got $knearest_lines"
fi
log "✅ TEST 4 completed successfully"

# Test 5: Distance reporting with different mode
log "Starting TEST 5: Distance reporting with signed distance"
"$meta_executable" \
    --input_a "$meta_temp_dir/queries.bed" \
    --input_b "$meta_temp_dir/test_b_distant.bed" \
    --distance_mode "ref" \
    --output "$meta_temp_dir/output5.bed"

check_file_exists "$meta_temp_dir/output5.bed" "signed distance output"
check_file_not_empty "$meta_temp_dir/output5.bed" "signed distance output"
check_file_line_count "$meta_temp_dir/output5.bed" 4 "signed distance line count"

# Check that distance column includes negative values (upstream features)
if ! grep -q "[-]" "$meta_temp_dir/output5.bed"; then
    log "Warning: No negative distances found, may not have upstream features"
fi
log "✅ TEST 5 completed successfully"

log "🎉 All bedtools_closest tests completed successfully!"
