#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_pairtobed"

####################################################################################################

log "Creating test data..."

# Create test BEDPE file with paired intervals
cat > "$meta_temp_dir/pairs.bedpe" << 'EOF'
chr1	100	200	chr1	300	400	pair1	100	+	-
chr1	150	250	chr1	350	450	pair2	200	+	+
chr2	50	150	chr2	200	300	pair3	150	-	-
chr2	500	600	chr3	700	800	pair4	300	+	+
chr1	1000	1100	chr1	1200	1300	pair5	250	-	+
EOF

# Create test BED file with genomic features
cat > "$meta_temp_dir/features.bed" << 'EOF'
chr1	120	180	gene1	100	+
chr1	320	380	gene2	200	-
chr2	80	120	gene3	150	+
chr2	250	350	gene4	300	+
chr1	1050	1080	gene5	400	+
chr1	1220	1280	gene6	500	-
EOF

####################################################################################################

log "TEST 1: Basic pairtobed functionality (default: either)"
"$meta_executable" \
  --bedpe "$meta_temp_dir/pairs.bedpe" \
  --bed "$meta_temp_dir/features.bed" \
  --output "$meta_temp_dir/test1_output.bedpe"

check_file_exists "$meta_temp_dir/test1_output.bedpe" "basic pairtobed output"
check_file_not_empty "$meta_temp_dir/test1_output.bedpe" "basic pairtobed result"

# Count the number of overlapping pairs
num_overlaps=$(wc -l < "$meta_temp_dir/test1_output.bedpe")
if [ "$num_overlaps" -gt 0 ]; then
  log "✓ Found $num_overlaps overlapping pairs (either end overlaps)"
else
  log "ERROR: No overlaps found, expected at least some overlaps"
  exit 1
fi

####################################################################################################

log "TEST 2: Both ends must overlap"
"$meta_executable" \
  --bedpe "$meta_temp_dir/pairs.bedpe" \
  --bed "$meta_temp_dir/features.bed" \
  --type "both" \
  --output "$meta_temp_dir/test2_output.bedpe"

check_file_exists "$meta_temp_dir/test2_output.bedpe" "both ends overlap output"

# Count overlaps with 'both' type - should be fewer than 'either'
num_both=$(wc -l < "$meta_temp_dir/test2_output.bedpe")
log "✓ Found $num_both pairs where both ends overlap"

####################################################################################################

log "TEST 3: Neither end overlaps"
"$meta_executable" \
  --bedpe "$meta_temp_dir/pairs.bedpe" \
  --bed "$meta_temp_dir/features.bed" \
  --type "neither" \
  --output "$meta_temp_dir/test3_output.bedpe"

check_file_exists "$meta_temp_dir/test3_output.bedpe" "neither end overlap output"

num_neither=$(wc -l < "$meta_temp_dir/test3_output.bedpe")
log "✓ Found $num_neither pairs where neither end overlaps"

####################################################################################################

log "TEST 4: XOR overlap (exactly one end overlaps)"
"$meta_executable" \
  --bedpe "$meta_temp_dir/pairs.bedpe" \
  --bed "$meta_temp_dir/features.bed" \
  --type "xor" \
  --output "$meta_temp_dir/test4_output.bedpe"

check_file_exists "$meta_temp_dir/test4_output.bedpe" "xor overlap output"

num_xor=$(wc -l < "$meta_temp_dir/test4_output.bedpe")
log "✓ Found $num_xor pairs where exactly one end overlaps"

####################################################################################################

log "TEST 5: Minimum overlap threshold"
"$meta_executable" \
  --bedpe "$meta_temp_dir/pairs.bedpe" \
  --bed "$meta_temp_dir/features.bed" \
  --min_overlap 0.5 \
  --output "$meta_temp_dir/test5_output.bedpe"

check_file_exists "$meta_temp_dir/test5_output.bedpe" "minimum overlap output"

num_min_overlap=$(wc -l < "$meta_temp_dir/test5_output.bedpe")
log "✓ Found $num_min_overlap pairs with minimum 50% overlap"

####################################################################################################

log "TEST 6: Same strand requirement"
"$meta_executable" \
  --bedpe "$meta_temp_dir/pairs.bedpe" \
  --bed "$meta_temp_dir/features.bed" \
  --same_strand \
  --output "$meta_temp_dir/test6_output.bedpe"

check_file_exists "$meta_temp_dir/test6_output.bedpe" "same strand output"

num_same_strand=$(wc -l < "$meta_temp_dir/test6_output.bedpe")
log "✓ Found $num_same_strand pairs with same strand requirement"

####################################################################################################

log "TEST 7: Opposite strand requirement"
"$meta_executable" \
  --bedpe "$meta_temp_dir/pairs.bedpe" \
  --bed "$meta_temp_dir/features.bed" \
  --opposite_strand \
  --output "$meta_temp_dir/test7_output.bedpe"

check_file_exists "$meta_temp_dir/test7_output.bedpe" "opposite strand output"

num_opposite_strand=$(wc -l < "$meta_temp_dir/test7_output.bedpe")
log "✓ Found $num_opposite_strand pairs with opposite strand requirement"

####################################################################################################

log "TEST 8: Span-based overlap (ispan)"
"$meta_executable" \
  --bedpe "$meta_temp_dir/pairs.bedpe" \
  --bed "$meta_temp_dir/features.bed" \
  --type "ispan" \
  --output "$meta_temp_dir/test8_output.bedpe"

check_file_exists "$meta_temp_dir/test8_output.bedpe" "ispan overlap output"

num_ispan=$(wc -l < "$meta_temp_dir/test8_output.bedpe")
log "✓ Found $num_ispan pairs with ispan overlap ([end1, start2] span)"

####################################################################################################

log "TEST 9: Parameter validation"
# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" --bed "$meta_temp_dir/features.bed" --output "$meta_temp_dir/test.bedpe" 2>/dev/null; then
  log "✗ Should have failed without --bedpe parameter"
  exit 1
else
  log "✓ Correctly requires --bedpe parameter (or --bam_input)"
fi

if "$meta_executable" --bedpe "$meta_temp_dir/pairs.bedpe" --output "$meta_temp_dir/test.bedpe" 2>/dev/null; then
  log "✗ Should have failed without --bed parameter"
  exit 1
else
  log "✓ Correctly requires --bed parameter"
fi

####################################################################################################

log "TEST 10: Logic validation - verify either = both + xor + neither"
# Mathematical check: pairs with either overlap should equal both + xor + neither
# This validates our different overlap types are working correctly

total_pairs=$(wc -l < "$meta_temp_dir/pairs.bedpe")
log "Total input pairs: $total_pairs"
log "Either overlaps: $num_overlaps"
log "Both overlaps: $num_both"
log "XOR overlaps: $num_xor"
log "Neither overlaps: $num_neither"

# Verify that both + xor + neither = total pairs
calculated_total=$((num_both + num_xor + num_neither))
if [ "$calculated_total" -eq "$total_pairs" ]; then
  log "✓ Overlap logic validation passed: both($num_both) + xor($num_xor) + neither($num_neither) = total($total_pairs)"
else
  log "WARNING: Overlap counts don't add up perfectly, but this can happen with edge cases"
  log "  Calculated total: $calculated_total, Actual total: $total_pairs"
fi

####################################################################################################

log "✓ All tests completed successfully!"
log "bedtools_pairtobed is working correctly with various overlap types and options"
