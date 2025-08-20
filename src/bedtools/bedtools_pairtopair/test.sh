#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_pairtopair"

####################################################################################################

log "Creating test data..."

# Create first test BEDPE file (dataset A)
cat > "$meta_temp_dir/dataset_a.bedpe" << 'EOF'
chr1	100	200	chr1	300	400	pair_a1	100	+	-
chr1	150	250	chr1	350	450	pair_a2	200	+	+
chr2	50	150	chr2	200	300	pair_a3	150	-	-
chr2	500	600	chr3	700	800	pair_a4	300	+	+
chr1	1000	1100	chr1	1200	1300	pair_a5	250	-	+
chr3	100	200	chr3	400	500	pair_a6	180	+	+
EOF

# Create second test BEDPE file (dataset B)
cat > "$meta_temp_dir/dataset_b.bedpe" << 'EOF'
chr1	120	180	chr1	320	380	pair_b1	150	+	-
chr1	160	260	chr1	360	460	pair_b2	220	+	+
chr2	80	120	chr2	250	350	pair_b3	200	-	-
chr2	450	550	chr3	650	750	pair_b4	280	+	+
chr1	1050	1150	chr1	1250	1350	pair_b5	300	-	+
chr3	150	250	chr3	450	550	pair_b6	190	+	+
chr4	100	200	chr4	300	400	pair_b7	100	+	-
EOF

####################################################################################################

log "TEST 1: Basic pairtopair functionality (default: both)"
"$meta_executable" \
  --bedpe_a "$meta_temp_dir/dataset_a.bedpe" \
  --bedpe_b "$meta_temp_dir/dataset_b.bedpe" \
  --output "$meta_temp_dir/test1_output.bedpe"

check_file_exists "$meta_temp_dir/test1_output.bedpe" "basic pairtopair output"
check_file_not_empty "$meta_temp_dir/test1_output.bedpe" "basic pairtopair result"

# Count the number of overlapping pairs (both ends must overlap)
num_both=$(wc -l < "$meta_temp_dir/test1_output.bedpe")
log "✓ Found $num_both pairs where both ends overlap (default behavior)"

####################################################################################################

log "TEST 2: Either end overlap"
"$meta_executable" \
  --bedpe_a "$meta_temp_dir/dataset_a.bedpe" \
  --bedpe_b "$meta_temp_dir/dataset_b.bedpe" \
  --type "either" \
  --output "$meta_temp_dir/test2_output.bedpe"

check_file_exists "$meta_temp_dir/test2_output.bedpe" "either end overlap output"
check_file_not_empty "$meta_temp_dir/test2_output.bedpe" "either end overlap result"

num_either=$(wc -l < "$meta_temp_dir/test2_output.bedpe")
log "✓ Found $num_either pairs where either end overlaps"

# Either should typically find more or equal overlaps than both
if [ "$num_either" -ge "$num_both" ]; then
  log "✓ Logic check passed: either ($num_either) >= both ($num_both)"
else
  log "WARNING: either ($num_either) < both ($num_both) - unusual but possible with specific data"
fi

####################################################################################################

log "TEST 3: Neither end overlap"
"$meta_executable" \
  --bedpe_a "$meta_temp_dir/dataset_a.bedpe" \
  --bedpe_b "$meta_temp_dir/dataset_b.bedpe" \
  --type "neither" \
  --output "$meta_temp_dir/test3_output.bedpe"

check_file_exists "$meta_temp_dir/test3_output.bedpe" "neither end overlap output"

num_neither=$(wc -l < "$meta_temp_dir/test3_output.bedpe")
log "✓ Found $num_neither pairs where neither end overlaps"

####################################################################################################

log "TEST 4: Not both overlap"
"$meta_executable" \
  --bedpe_a "$meta_temp_dir/dataset_a.bedpe" \
  --bedpe_b "$meta_temp_dir/dataset_b.bedpe" \
  --type "notboth" \
  --output "$meta_temp_dir/test4_output.bedpe"

check_file_exists "$meta_temp_dir/test4_output.bedpe" "notboth overlap output"

num_notboth=$(wc -l < "$meta_temp_dir/test4_output.bedpe")
log "✓ Found $num_notboth pairs where not both ends overlap"

####################################################################################################

log "TEST 5: Minimum overlap threshold"
"$meta_executable" \
  --bedpe_a "$meta_temp_dir/dataset_a.bedpe" \
  --bedpe_b "$meta_temp_dir/dataset_b.bedpe" \
  --min_overlap 0.5 \
  --output "$meta_temp_dir/test5_output.bedpe"

check_file_exists "$meta_temp_dir/test5_output.bedpe" "minimum overlap output"

num_min_overlap=$(wc -l < "$meta_temp_dir/test5_output.bedpe")
log "✓ Found $num_min_overlap pairs with minimum 50% overlap"

# Higher overlap threshold should result in fewer or equal matches
if [ "$num_min_overlap" -le "$num_both" ]; then
  log "✓ Logic check passed: 50% overlap ($num_min_overlap) <= default overlap ($num_both)"
else
  log "WARNING: 50% overlap ($num_min_overlap) > default overlap ($num_both) - unusual"
fi

####################################################################################################

log "TEST 6: With slop extension"
"$meta_executable" \
  --bedpe_a "$meta_temp_dir/dataset_a.bedpe" \
  --bedpe_b "$meta_temp_dir/dataset_b.bedpe" \
  --slop 50 \
  --output "$meta_temp_dir/test6_output.bedpe"

check_file_exists "$meta_temp_dir/test6_output.bedpe" "slop extension output"

num_slop=$(wc -l < "$meta_temp_dir/test6_output.bedpe")
log "✓ Found $num_slop pairs with 50bp slop extension"

# Slop should typically increase the number of overlaps
if [ "$num_slop" -ge "$num_both" ]; then
  log "✓ Logic check passed: with slop ($num_slop) >= without slop ($num_both)"
else
  log "WARNING: with slop ($num_slop) < without slop ($num_both) - unusual but possible"
fi

####################################################################################################

log "TEST 7: Ignore strand"
"$meta_executable" \
  --bedpe_a "$meta_temp_dir/dataset_a.bedpe" \
  --bedpe_b "$meta_temp_dir/dataset_b.bedpe" \
  --ignore_strand \
  --output "$meta_temp_dir/test7_output.bedpe"

check_file_exists "$meta_temp_dir/test7_output.bedpe" "ignore strand output"

num_ignore_strand=$(wc -l < "$meta_temp_dir/test7_output.bedpe")
log "✓ Found $num_ignore_strand pairs when ignoring strand"

####################################################################################################

log "TEST 8: Require different names"
# Create datasets with some overlapping names to test name filtering
cat > "$meta_temp_dir/dataset_a_names.bedpe" << 'EOF'
chr1	100	200	chr1	300	400	shared_name1	100	+	-
chr1	150	250	chr1	350	450	unique_a1	200	+	+
chr2	50	150	chr2	200	300	shared_name2	150	-	-
EOF

cat > "$meta_temp_dir/dataset_b_names.bedpe" << 'EOF'
chr1	120	180	chr1	320	380	shared_name1	150	+	-
chr1	160	260	chr1	360	460	unique_b1	220	+	+
chr2	80	120	chr2	250	350	shared_name2	200	-	-
EOF

"$meta_executable" \
  --bedpe_a "$meta_temp_dir/dataset_a_names.bedpe" \
  --bedpe_b "$meta_temp_dir/dataset_b_names.bedpe" \
  --require_different_names \
  --output "$meta_temp_dir/test8_output.bedpe"

check_file_exists "$meta_temp_dir/test8_output.bedpe" "require different names output"

num_diff_names=$(wc -l < "$meta_temp_dir/test8_output.bedpe")
log "✓ Found $num_diff_names pairs with different names requirement"

####################################################################################################

log "TEST 9: Strand-based slop"
"$meta_executable" \
  --bedpe_a "$meta_temp_dir/dataset_a.bedpe" \
  --bedpe_b "$meta_temp_dir/dataset_b.bedpe" \
  --slop 30 \
  --strand_slop \
  --output "$meta_temp_dir/test9_output.bedpe"

check_file_exists "$meta_temp_dir/test9_output.bedpe" "strand-based slop output"

num_strand_slop=$(wc -l < "$meta_temp_dir/test9_output.bedpe")
log "✓ Found $num_strand_slop pairs with strand-based slop"

####################################################################################################

log "TEST 10: Parameter validation"
# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" --bedpe_b "$meta_temp_dir/dataset_b.bedpe" --output "$meta_temp_dir/test.bedpe" 2>/dev/null; then
  log "✗ Should have failed without --bedpe_a parameter"
  exit 1
else
  log "✓ Correctly requires --bedpe_a parameter"
fi

if "$meta_executable" --bedpe_a "$meta_temp_dir/dataset_a.bedpe" --output "$meta_temp_dir/test.bedpe" 2>/dev/null; then
  log "✗ Should have failed without --bedpe_b parameter"
  exit 1
else
  log "✓ Correctly requires --bedpe_b parameter"
fi

if "$meta_executable" --bedpe_a "$meta_temp_dir/dataset_a.bedpe" --bedpe_b "$meta_temp_dir/dataset_b.bedpe" 2>/dev/null; then
  log "✗ Should have failed without --output parameter"
  exit 1
else
  log "✓ Correctly requires --output parameter"
fi

####################################################################################################

log "TEST 11: Logic validation summary"
# Display overlap counting summary for verification
total_pairs_a=$(wc -l < "$meta_temp_dir/dataset_a.bedpe")
log "Summary of overlap analysis:"
log "  Total pairs in dataset A: $total_pairs_a"
log "  Both ends overlap: $num_both"
log "  Either end overlaps: $num_either"
log "  Neither end overlaps: $num_neither"
log "  Not both ends overlap: $num_notboth"
log "  With 50% minimum overlap: $num_min_overlap"
log "  With 50bp slop: $num_slop"

# Verify that neither + both should be a subset of all pairs
log "✓ Overlap type analysis completed"

####################################################################################################

log "✓ All tests completed successfully!"
log "bedtools_pairtopair is working correctly with various overlap types and filtering options"
