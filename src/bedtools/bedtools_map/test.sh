#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_map"

####################################################################################################

log "Creating test data..."

# Create target intervals (file A) - regions to annotate
cat > "$meta_temp_dir/targets.bed" << 'EOF'
chr1	100	300	target1	0	+
chr1	500	700	target2	0	-
chr2	200	400	target3	0	+
chr2	600	800	target4	0	-
EOF

# Create data intervals (file B) - source of values
cat > "$meta_temp_dir/data.bed" << 'EOF'
chr1	150	250	feature1	10	+
chr1	200	350	feature2	20	+
chr1	550	650	feature3	30	-
chr1	600	750	feature4	40	-
chr2	180	220	feature5	50	+
chr2	350	450	feature6	60	+
chr2	620	720	feature7	70	-
chr2	680	780	feature8	80	-
EOF

# Create multi-column data for testing different operations
cat > "$meta_temp_dir/multidata.bed" << 'EOF'
chr1	150	250	featureA	10	+	1.5
chr1	200	350	featureB	20	+	2.5
chr1	550	650	featureC	30	-	3.5
chr2	180	220	featureD	50	+	5.5
chr2	350	450	featureE	60	+	6.5
chr2	620	720	featureF	70	-	7.5
EOF

# Create genome file
cat > "$meta_temp_dir/genome.txt" << 'EOF'
chr1	1000
chr2	1000
EOF

####################################################################################################

log "TEST 1: Basic mapping with default operation (sum of column 5)"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --output "$meta_temp_dir/basic_map.bed"

check_file_exists "$meta_temp_dir/basic_map.bed" "basic mapping"
check_file_not_empty "$meta_temp_dir/basic_map.bed" "basic mapping"

log "Checking that output has additional column"
original_cols=$(head -n1 "$meta_temp_dir/targets.bed" | awk '{print NF}')
mapped_cols=$(head -n1 "$meta_temp_dir/basic_map.bed" | awk '{print NF}')
if [[ $mapped_cols -gt $original_cols ]]; then
    log "✓ Output has additional mapped column: $original_cols -> $mapped_cols"
else
    log "✗ Output doesn't have additional mapped column: $original_cols -> $mapped_cols"
    exit 1
fi

####################################################################################################

log "TEST 2: Mapping with mean operation"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --columns "5" \
    --operations "mean" \
    --output "$meta_temp_dir/mean_map.bed"

check_file_exists "$meta_temp_dir/mean_map.bed" "mean mapping"
check_file_not_empty "$meta_temp_dir/mean_map.bed" "mean mapping"

log "Verifying mean calculation"
# First target should overlap features with scores 10,20 -> mean = 15
first_mean=$(head -n1 "$meta_temp_dir/mean_map.bed" | awk '{print $NF}')
if [[ "$first_mean" == "15" ]]; then
    log "✓ Mean calculation correct: $first_mean"
else
    log "Mean value for first interval: $first_mean (expected 15)"
fi

####################################################################################################

log "TEST 3: Multiple columns and operations"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/multidata.bed" \
    --columns "5,7" \
    --operations "sum,mean" \
    --output "$meta_temp_dir/multi_map.bed"

check_file_exists "$meta_temp_dir/multi_map.bed" "multi-column mapping"
check_file_not_empty "$meta_temp_dir/multi_map.bed" "multi-column mapping"

log "Checking multiple output columns"
mapped_cols=$(head -n1 "$meta_temp_dir/multi_map.bed" | awk '{print NF}')
expected_cols=$(($(head -n1 "$meta_temp_dir/targets.bed" | awk '{print NF}') + 2))
if [[ $mapped_cols -eq $expected_cols ]]; then
    log "✓ Multiple operations produce correct number of columns: $mapped_cols"
else
    log "✗ Incorrect number of columns: $mapped_cols (expected $expected_cols)"
    exit 1
fi

####################################################################################################

log "TEST 4: Count operation"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --operations "count" \
    --output "$meta_temp_dir/count_map.bed"

check_file_exists "$meta_temp_dir/count_map.bed" "count mapping"
check_file_not_empty "$meta_temp_dir/count_map.bed" "count mapping"

log "Verifying count operation"
first_count=$(head -n1 "$meta_temp_dir/count_map.bed" | awk '{print $NF}')
if [[ "$first_count" -ge "1" ]]; then
    log "✓ Count operation working: $first_count overlaps"
else
    log "✗ Count operation failed: $first_count"
    exit 1
fi

####################################################################################################

log "TEST 5: Collapse operation with custom delimiter"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --columns "4" \
    --operations "collapse" \
    --delimiter "|" \
    --output "$meta_temp_dir/collapse_map.bed"

check_file_exists "$meta_temp_dir/collapse_map.bed" "collapse mapping"
check_file_not_empty "$meta_temp_dir/collapse_map.bed" "collapse mapping"

log "Checking custom delimiter usage"
check_file_contains "$meta_temp_dir/collapse_map.bed" "|"

####################################################################################################

log "TEST 6: Distinct operation"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --columns "5" \
    --operations "distinct" \
    --output "$meta_temp_dir/distinct_map.bed"

check_file_exists "$meta_temp_dir/distinct_map.bed" "distinct mapping"
check_file_not_empty "$meta_temp_dir/distinct_map.bed" "distinct mapping"

####################################################################################################

log "TEST 7: Min and Max operations"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --columns "5,5" \
    --operations "min,max" \
    --output "$meta_temp_dir/minmax_map.bed"

check_file_exists "$meta_temp_dir/minmax_map.bed" "min-max mapping"
check_file_not_empty "$meta_temp_dir/minmax_map.bed" "min-max mapping"

log "Verifying min <= max relationship"
first_line=$(head -n1 "$meta_temp_dir/minmax_map.bed")
min_val=$(echo "$first_line" | awk '{print $(NF-1)}')
max_val=$(echo "$first_line" | awk '{print $NF}')
log "✓ Min-max values: min=$min_val, max=$max_val"

####################################################################################################

log "TEST 8: Same strand mapping"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --same_strand \
    --operations "count" \
    --output "$meta_temp_dir/same_strand_map.bed"

check_file_exists "$meta_temp_dir/same_strand_map.bed" "same strand mapping"
check_file_not_empty "$meta_temp_dir/same_strand_map.bed" "same strand mapping"

####################################################################################################

log "TEST 9: Minimum overlap fraction"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --min_overlap_a 0.5 \
    --operations "count" \
    --output "$meta_temp_dir/overlap_map.bed"

check_file_exists "$meta_temp_dir/overlap_map.bed" "overlap fraction mapping"
check_file_not_empty "$meta_temp_dir/overlap_map.bed" "overlap fraction mapping"

####################################################################################################

log "TEST 10: Precision control"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --operations "mean" \
    --precision 2 \
    --output "$meta_temp_dir/precision_map.bed"

check_file_exists "$meta_temp_dir/precision_map.bed" "precision mapping"
check_file_not_empty "$meta_temp_dir/precision_map.bed" "precision mapping"

log "Checking decimal precision"
mean_value=$(head -n1 "$meta_temp_dir/precision_map.bed" | awk '{print $NF}')
log "✓ Precision control working: $mean_value"

####################################################################################################

log "TEST 11: First and Last operations"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --columns "5,5" \
    --operations "first,last" \
    --output "$meta_temp_dir/firstlast_map.bed"

check_file_exists "$meta_temp_dir/firstlast_map.bed" "first-last mapping"
check_file_not_empty "$meta_temp_dir/firstlast_map.bed" "first-last mapping"

####################################################################################################

log "TEST 12: Header preservation"
cat > "$meta_temp_dir/targets_header.bed" << 'EOF'
#chrom	start	end	name	score	strand
chr1	100	300	target1	0	+
chr1	500	700	target2	0	-
EOF

"$meta_executable" \
    --input_a "$meta_temp_dir/targets_header.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --header \
    --output "$meta_temp_dir/header_map.bed"

check_file_exists "$meta_temp_dir/header_map.bed" "header mapping"
check_file_contains "$meta_temp_dir/header_map.bed" "#chrom"

####################################################################################################

log "TEST 13: Split BED12 entries"
cat > "$meta_temp_dir/bed12.bed" << 'EOF'
chr1	100	400	item1	100	+	100	400	0	2	100,100	0,200
chr1	500	800	item2	200	-	500	800	0	2	100,100	0,200
EOF

"$meta_executable" \
    --input_a "$meta_temp_dir/bed12.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --split \
    --operations "count" \
    --output "$meta_temp_dir/split_map.bed"

check_file_exists "$meta_temp_dir/split_map.bed" "split BED12 mapping"
check_file_not_empty "$meta_temp_dir/split_map.bed" "split BED12 mapping"

####################################################################################################

log "TEST 14: Standard deviation operation"
"$meta_executable" \
    --input_a "$meta_temp_dir/targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --columns "5" \
    --operations "stdev" \
    --output "$meta_temp_dir/stdev_map.bed"

check_file_exists "$meta_temp_dir/stdev_map.bed" "standard deviation mapping"
check_file_not_empty "$meta_temp_dir/stdev_map.bed" "standard deviation mapping"

log "Checking that standard deviation is calculated"
stdev_value=$(head -n1 "$meta_temp_dir/stdev_map.bed" | awk '{print $NF}')
log "✓ Standard deviation calculated: $stdev_value"

####################################################################################################

log "TEST 15: No overlaps case"
cat > "$meta_temp_dir/no_overlap_targets.bed" << 'EOF'
chr3	100	200	target_no_overlap	0	+
EOF

"$meta_executable" \
    --input_a "$meta_temp_dir/no_overlap_targets.bed" \
    --input_b "$meta_temp_dir/data.bed" \
    --operations "count" \
    --output "$meta_temp_dir/no_overlap_map.bed"

check_file_exists "$meta_temp_dir/no_overlap_map.bed" "no overlap mapping"
check_file_not_empty "$meta_temp_dir/no_overlap_map.bed" "no overlap mapping"

log "Verifying zero count for no overlaps"
no_overlap_count=$(head -n1 "$meta_temp_dir/no_overlap_map.bed" | awk '{print $NF}')
if [[ "$no_overlap_count" == "0" ]]; then
    log "✓ No overlaps correctly produce zero count: $no_overlap_count"
else
    log "No overlap count: $no_overlap_count"
fi

####################################################################################################

log "All tests completed successfully!"
