#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_overlap"

####################################################################################################

log "Creating test data..."

# Create test input file similar to bedtools window output
# This simulates the output from: bedtools window -a A.bed -b B.bed -w 10
cat > "$meta_temp_dir/windowed_features.bed" << 'EOF'
chr1	10	20	A	chr1	15	25	B
chr1	10	20	C	chr1	25	35	D
chr1	100	200	E	chr1	150	250	F
chr2	50	100	G	chr2	80	120	H
chr2	300	400	I	chr2	450	500	J
EOF

####################################################################################################

log "TEST 1: Basic overlap functionality"
"$meta_executable" \
  --input "$meta_temp_dir/windowed_features.bed" \
  --cols "2,3,6,7" \
  --output "$meta_temp_dir/test1_output.bed"

check_file_exists "$meta_temp_dir/test1_output.bed" "basic overlap output"
check_file_not_empty "$meta_temp_dir/test1_output.bed" "basic overlap result"

# Verify that output has the correct number of columns (original + 1 overlap column)
num_columns=$(head -1 "$meta_temp_dir/test1_output.bed" | awk '{print NF}')
if [ "$num_columns" -ne 9 ]; then  # 8 original columns + 1 overlap column
  log "ERROR: Output should have 9 columns (8 original + 1 overlap), found $num_columns"
  head -3 "$meta_temp_dir/test1_output.bed"
  exit 1
fi
log "✓ Output has correct number of columns ($num_columns)"

# Check that overlap values are computed correctly
# Line 1: chr1 10-20 vs chr1 15-25 should have overlap of 5
# Line 2: chr1 10-20 vs chr1 25-35 should have distance of -5 (negative)
overlap1=$(head -1 "$meta_temp_dir/test1_output.bed" | awk '{print $9}')
overlap2=$(sed -n '2p' "$meta_temp_dir/test1_output.bed" | awk '{print $9}')

if [ "$overlap1" = "5" ]; then
  log "✓ First overlap calculation correct: $overlap1"
else
  log "ERROR: Expected overlap of 5 for first line, got: $overlap1"
  exit 1
fi

if [ "$overlap2" = "-5" ]; then
  log "✓ Second distance calculation correct: $overlap2"
else
  log "ERROR: Expected distance of -5 for second line, got: $overlap2"
  exit 1
fi

####################################################################################################

log "TEST 2: Different column specification"
# Test with different column positions
cat > "$meta_temp_dir/custom_format.bed" << 'EOF'
feature1	chr1	100	150	feature2	chr1	120	180
feature3	chr1	200	250	feature4	chr1	300	350
EOF

"$meta_executable" \
  --input "$meta_temp_dir/custom_format.bed" \
  --cols "3,4,7,8" \
  --output "$meta_temp_dir/test2_output.bed"

check_file_exists "$meta_temp_dir/test2_output.bed" "custom columns output"
check_file_not_empty "$meta_temp_dir/test2_output.bed" "custom columns result"

# Check that the first line has correct overlap (100-150 vs 120-180 = 30 overlap)
overlap_custom=$(head -1 "$meta_temp_dir/test2_output.bed" | awk '{print $NF}')
if [ "$overlap_custom" = "30" ]; then
  log "✓ Custom column overlap calculation correct: $overlap_custom"
else
  log "ERROR: Expected overlap of 30 for custom columns, got: $overlap_custom"
  exit 1
fi

####################################################################################################

log "TEST 3: Multiple overlap scenarios"
# Test various overlap and distance scenarios
cat > "$meta_temp_dir/multiple_scenarios.bed" << 'EOF'
chr1	0	10	A	chr1	5	15	B
chr1	20	30	C	chr1	35	45	D  
chr1	50	100	E	chr1	40	60	F
chr1	200	300	G	chr1	200	300	H
EOF

"$meta_executable" \
  --input "$meta_temp_dir/multiple_scenarios.bed" \
  --cols "2,3,6,7" \
  --output "$meta_temp_dir/test3_output.bed"

check_file_exists "$meta_temp_dir/test3_output.bed" "multiple scenarios output"

# Verify various calculations:
# Line 1: 0-10 vs 5-15 = 5 overlap
# Line 2: 20-30 vs 35-45 = -5 distance  
# Line 3: 50-100 vs 40-60 = 10 overlap
# Line 4: 200-300 vs 200-300 = 100 overlap (identical)

overlaps=($(awk '{print $NF}' "$meta_temp_dir/test3_output.bed"))
expected=(5 -5 10 100)

for i in {0..3}; do
  if [ "${overlaps[i]}" = "${expected[i]}" ]; then
      log "✓ Scenario $((i+1)) overlap correct: ${overlaps[i]}"
  else
      log "ERROR: Scenario $((i+1)) expected ${expected[i]}, got: ${overlaps[i]}"
      exit 1
  fi
done

####################################################################################################

log "TEST 4: Parameter validation"
# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" --input "$meta_temp_dir/windowed_features.bed" --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "✗ Should have failed without --cols parameter"
  exit 1
else
  log "✓ Correctly requires --cols parameter"
fi

if "$meta_executable" --cols "2,3,6,7" --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "✗ Should have failed without --input parameter"
  exit 1
else
  log "✓ Correctly requires --input parameter"
fi

####################################################################################################

log "✓ All tests completed successfully!"
log "bedtools_overlap is working correctly with proper overlap and distance calculations"
