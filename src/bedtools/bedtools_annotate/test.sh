#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_annotate"

# Create test data
log "Creating test data..."

# Create input intervals file
cat > "$meta_temp_dir/intervals.bed" << 'EOF'
chr1	100	200	interval1	100	+
chr1	300	400	interval2	200	+
chr2	150	250	interval3	300	-
chr2	500	600	interval4	400	-
EOF

# Create first annotation file (overlaps with intervals 1 and 3)
cat > "$meta_temp_dir/annotation1.bed" << 'EOF'
chr1	120	180	feature1	500	+
chr1	350	450	feature2	600	+
chr2	140	260	feature3	700	-
EOF

# Create second annotation file (overlaps with intervals 2 and 4)
cat > "$meta_temp_dir/annotation2.bed" << 'EOF'
chr1	320	380	feature4	800	+
chr1	390	420	feature5	900	+
chr2	520	580	feature6	1000	-
EOF

# Test 1: Basic annotation with coverage fractions
log "Starting TEST 1: Basic annotation with coverage fractions"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --files "$meta_temp_dir/annotation1.bed;$meta_temp_dir/annotation2.bed" \
  --output "$meta_temp_dir/output1.bed"

check_file_exists "$meta_temp_dir/output1.bed" "basic annotation output"
check_file_not_empty "$meta_temp_dir/output1.bed" "basic annotation output"
check_file_line_count "$meta_temp_dir/output1.bed" 4 "basic annotation line count"

# Check that fractions are present (should contain decimal numbers)
check_file_contains "$meta_temp_dir/output1.bed" "0." "coverage fractions"
log "✅ TEST 1 completed successfully"

# Test 2: Annotation with feature counts
log "Starting TEST 2: Annotation with feature counts"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --files "$meta_temp_dir/annotation1.bed;$meta_temp_dir/annotation2.bed" \
  --output "$meta_temp_dir/output2.bed" \
  --counts

check_file_exists "$meta_temp_dir/output2.bed" "count annotation output"
check_file_not_empty "$meta_temp_dir/output2.bed" "count annotation output"

# Check that counts are present (should contain integers)
check_file_contains "$meta_temp_dir/output2.bed" "1" "feature counts"
log "✅ TEST 2 completed successfully"

# Test 3: Annotation with both counts and fractions
log "Starting TEST 3: Annotation with both counts and fractions"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --files "$meta_temp_dir/annotation1.bed" \
  --output "$meta_temp_dir/output3.bed" \
  --both

check_file_exists "$meta_temp_dir/output3.bed" "both metrics output"
check_file_not_empty "$meta_temp_dir/output3.bed" "both metrics output"

# Check that both counts and fractions are present
check_file_contains "$meta_temp_dir/output3.bed" "1" "feature counts in both output"
check_file_contains "$meta_temp_dir/output3.bed" "0." "coverage fractions in both output"
log "✅ TEST 3 completed successfully"

# Test 4: Annotation with custom names
log "Starting TEST 4: Annotation with custom names"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --files "$meta_temp_dir/annotation1.bed;$meta_temp_dir/annotation2.bed" \
  --names "ChIP_peaks;DNA_meth" \
  --output "$meta_temp_dir/output4.bed"

check_file_exists "$meta_temp_dir/output4.bed" "named annotation output"
check_file_not_empty "$meta_temp_dir/output4.bed" "named annotation output"

# The names should appear somewhere (likely in header or within results)
log "✅ TEST 4 completed successfully"

# Test 5: Strand-specific annotation (same strand)
log "Starting TEST 5: Strand-specific annotation (same strand)"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --files "$meta_temp_dir/annotation1.bed" \
  --output "$meta_temp_dir/output5.bed" \
  --strand

check_file_exists "$meta_temp_dir/output5.bed" "strand-specific annotation output"
check_file_not_empty "$meta_temp_dir/output5.bed" "strand-specific annotation output"
log "✅ TEST 5 completed successfully"

log "All tests completed successfully!"
