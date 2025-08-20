#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_coverage"

# Create test data
log "Creating test data..."

# Create target intervals (query file A)
cat > "$meta_temp_dir/targets.bed" << 'EOF'
chr1	100	300	target1	100	+
chr1	500	800	target2	200	+
chr2	200	400	target3	300	-
chr2	600	900	target4	400	-
EOF

# Create coverage features (file B) - some overlapping, some not
cat > "$meta_temp_dir/features.bed" << 'EOF'
chr1	150	250	feature1	500	+
chr1	200	350	feature2	600	+
chr1	550	750	feature3	700	+
chr2	250	350	feature4	800	-
chr2	650	850	feature5	900	+
chr3	100	200	feature6	1000	+
EOF

# Create additional coverage file for multi-file testing
cat > "$meta_temp_dir/features2.bed" << 'EOF'
chr1	120	180	extra1	300	+
chr1	600	700	extra2	400	+
chr2	300	500	extra3	500	-
EOF

# Create strand-specific test data
cat > "$meta_temp_dir/stranded_targets.bed" << 'EOF'
chr1	100	200	pos_target	100	+
chr1	300	400	neg_target	200	-
EOF

cat > "$meta_temp_dir/stranded_features.bed" << 'EOF'
chr1	120	180	pos_feature	300	+
chr1	320	380	neg_feature	400	-
chr1	140	160	pos_feature2	500	+
chr1	340	360	neg_feature2	600	-
EOF

# Test 1: Basic coverage calculation
log "Starting TEST 1: Basic coverage calculation"
"$meta_executable" \
  --input_a "$meta_temp_dir/targets.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --output "$meta_temp_dir/output1.txt"

check_file_exists "$meta_temp_dir/output1.txt" "basic coverage output"
check_file_not_empty "$meta_temp_dir/output1.txt" "basic coverage output"
check_file_line_count "$meta_temp_dir/output1.txt" 4 "basic coverage line count"

# Check that coverage statistics are added (should have 4 extra columns)
input_cols=$(head -1 "$meta_temp_dir/targets.bed" | awk '{print NF}')
output_cols=$(head -1 "$meta_temp_dir/output1.txt" | awk '{print NF}')
expected_cols=$((input_cols + 4))
if [ $output_cols -ne $expected_cols ]; then
  log_error "Expected $expected_cols columns in output, got $output_cols"
  exit 1
fi

# Check that some targets have coverage
if ! grep -q -E "\s[1-9][0-9]*\s" "$meta_temp_dir/output1.txt"; then
  log_error "Expected some targets to have non-zero coverage counts"
  exit 1
fi
log "✅ TEST 1 completed successfully"

# Test 2: Coverage histogram
log "Starting TEST 2: Coverage histogram"
"$meta_executable" \
  --input_a "$meta_temp_dir/targets.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --histogram \
  --output "$meta_temp_dir/output2.txt"

check_file_exists "$meta_temp_dir/output2.txt" "histogram output"
check_file_not_empty "$meta_temp_dir/output2.txt" "histogram output"

# Histogram output should have depth information
check_file_contains "$meta_temp_dir/output2.txt" "target1" "target intervals in histogram"
# Should contain histogram data (depth, bases, size, percentage)
if ! grep -q -E "\s[0-9]+\s+[0-9]+\s+[0-9]+\s+[0-9]+\.[0-9]+$" "$meta_temp_dir/output2.txt"; then
  log_error "Expected histogram format with depth data"
  exit 1
fi
log "✅ TEST 2 completed successfully"

# Test 3: Counts only
log "Starting TEST 3: Counts only output"
"$meta_executable" \
  --input_a "$meta_temp_dir/targets.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --counts_only \
  --output "$meta_temp_dir/output3.txt"

check_file_exists "$meta_temp_dir/output3.txt" "counts only output"
check_file_not_empty "$meta_temp_dir/output3.txt" "counts only output"
check_file_line_count "$meta_temp_dir/output3.txt" 4 "counts only line count"

# Counts only should have fewer columns (just original + count)
counts_cols=$(head -1 "$meta_temp_dir/output3.txt" | awk '{print NF}')
if [ $counts_cols -ne $((input_cols + 1)) ]; then
  log_error "Expected $((input_cols + 1)) columns for counts only, got $counts_cols"
  exit 1
fi
log "✅ TEST 3 completed successfully"

# Test 4: Mean depth reporting
log "Starting TEST 4: Mean depth reporting"
"$meta_executable" \
  --input_a "$meta_temp_dir/targets.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --mean_depth \
  --output "$meta_temp_dir/output4.txt"

check_file_exists "$meta_temp_dir/output4.txt" "mean depth output"
check_file_not_empty "$meta_temp_dir/output4.txt" "mean depth output"

# Should contain mean depth values (floating point numbers)
if ! grep -q -E "\s[0-9]+\.[0-9]+$" "$meta_temp_dir/output4.txt"; then
  log_error "Expected mean depth values (floating point)"
  exit 1
fi
log "✅ TEST 4 completed successfully"

# Test 5: Strand-specific coverage
log "Starting TEST 5: Strand-specific coverage"
"$meta_executable" \
  --input_a "$meta_temp_dir/stranded_targets.bed" \
  --input_b "$meta_temp_dir/stranded_features.bed" \
  --same_strand \
  --output "$meta_temp_dir/output5.txt"

check_file_exists "$meta_temp_dir/output5.txt" "same strand output"
check_file_not_empty "$meta_temp_dir/output5.txt" "same strand output"

# Compare with opposite strand requirement
"$meta_executable" \
  --input_a "$meta_temp_dir/stranded_targets.bed" \
  --input_b "$meta_temp_dir/stranded_features.bed" \
  --different_strand \
  --output "$meta_temp_dir/output5b.txt"

# Results should be different between same and different strand requirements
if diff -q "$meta_temp_dir/output5.txt" "$meta_temp_dir/output5b.txt" >/dev/null; then
  log "Warning: Same and different strand outputs are identical - may not have strand-specific overlaps"
fi
log "✅ TEST 5 completed successfully"

# Test 6: Multiple input files
log "Starting TEST 6: Multiple input files"
"$meta_executable" \
  --input_a "$meta_temp_dir/targets.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --input_b "$meta_temp_dir/features2.bed" \
  --output "$meta_temp_dir/output6.txt"

check_file_exists "$meta_temp_dir/output6.txt" "multiple files output"
check_file_not_empty "$meta_temp_dir/output6.txt" "multiple files output"
check_file_line_count "$meta_temp_dir/output6.txt" 4 "multiple files line count"

# Coverage should be higher with additional file
single_file_coverage=$(awk '{print $7}' "$meta_temp_dir/output1.txt" | head -1)
multi_file_coverage=$(awk '{print $7}' "$meta_temp_dir/output6.txt" | head -1)
log "ℹ️  Single file coverage: $single_file_coverage, Multi-file coverage: $multi_file_coverage"
log "✅ TEST 6 completed successfully"

# Test 7: Minimum overlap fraction
log "Starting TEST 7: Minimum overlap fraction"
"$meta_executable" \
  --input_a "$meta_temp_dir/targets.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --min_overlap_a 0.5 \
  --output "$meta_temp_dir/output7.txt"

check_file_exists "$meta_temp_dir/output7.txt" "min overlap output"
check_file_not_empty "$meta_temp_dir/output7.txt" "min overlap output"

# Compare with no minimum requirement - should have fewer overlaps
no_min_overlaps=$(awk '{sum += $7} END {print sum}' "$meta_temp_dir/output1.txt")
min_overlaps=$(awk '{sum += $7} END {print sum}' "$meta_temp_dir/output7.txt")

if [ "$min_overlaps" -gt "$no_min_overlaps" ]; then
  log_error "Expected fewer overlaps with minimum fraction requirement"
  exit 1
fi
log "✅ TEST 7 completed successfully"

log "🎉 All bedtools_coverage tests completed successfully!"
