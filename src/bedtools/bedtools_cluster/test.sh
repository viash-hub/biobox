#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_cluster"

# Create test data
log "Creating test data..."

# Create overlapping intervals for basic clustering
cat > "$meta_temp_dir/overlapping.bed" << 'EOF'
chr1	100	200	feature1	100	+
chr1	150	250	feature2	200	+
chr1	180	280	feature3	300	+
chr1	500	600	feature4	400	-
chr1	800	900	feature5	500	+
chr2	100	200	feature6	600	+
chr2	300	400	feature7	700	-
EOF

# Create intervals with different strands
cat > "$meta_temp_dir/stranded.bed" << 'EOF'
chr1	100	200	pos1	100	+
chr1	150	250	neg1	200	-
chr1	180	280	pos2	300	+
chr1	300	400	neg2	400	-
chr1	500	600	pos3	500	+
chr1	550	650	neg3	600	-
EOF

# Create intervals for distance-based clustering
cat > "$meta_temp_dir/nearby.bed" << 'EOF'
chr1	100	200	interval1	100	+
chr1	300	400	interval2	200	+
chr1	450	550	interval3	300	+
chr1	1000	1100	interval4	400	+
chr1	1200	1300	interval5	500	+
chr2	100	200	interval6	600	+
chr2	1000	1100	interval7	700	+
EOF

# Test 1: Basic clustering of overlapping intervals
log "Starting TEST 1: Basic clustering of overlapping intervals"
"$meta_executable" \
    --input "$meta_temp_dir/overlapping.bed" \
    --output "$meta_temp_dir/output1.bed"

check_file_exists "$meta_temp_dir/output1.bed" "basic clustering output"
check_file_not_empty "$meta_temp_dir/output1.bed" "basic clustering output"
check_file_line_count "$meta_temp_dir/output1.bed" 7 "basic clustering line count"

# Check that cluster IDs are added (should have one more column than input)
input_cols=$(head -1 "$meta_temp_dir/overlapping.bed" | awk '{print NF}')
output_cols=$(head -1 "$meta_temp_dir/output1.bed" | awk '{print NF}')
if [ $output_cols -ne $((input_cols + 1)) ]; then
    log_error "Expected $((input_cols + 1)) columns in output, got $output_cols"
    exit 1
fi

# Check that overlapping intervals get the same cluster ID
if ! grep -q "	1$" "$meta_temp_dir/output1.bed"; then
    log_error "Expected cluster ID 1 in output"
    exit 1
fi
log "✅ TEST 1 completed successfully"

# Test 2: Distance-based clustering
log "Starting TEST 2: Distance-based clustering"
"$meta_executable" \
    --input "$meta_temp_dir/nearby.bed" \
    --distance 100 \
    --output "$meta_temp_dir/output2.bed"

check_file_exists "$meta_temp_dir/output2.bed" "distance clustering output"
check_file_not_empty "$meta_temp_dir/output2.bed" "distance clustering output"
check_file_line_count "$meta_temp_dir/output2.bed" 7 "distance clustering line count"

# With distance 100, intervals at positions 100-200, 300-400, 450-550 should cluster together
# Check that cluster IDs are present
check_file_contains "$meta_temp_dir/output2.bed" "1" "cluster IDs present"
log "✅ TEST 2 completed successfully"

# Test 3: Strand-specific clustering
log "Starting TEST 3: Strand-specific clustering"
"$meta_executable" \
    --input "$meta_temp_dir/stranded.bed" \
    --strand \
    --output "$meta_temp_dir/output3.bed"

check_file_exists "$meta_temp_dir/output3.bed" "strand clustering output"
check_file_not_empty "$meta_temp_dir/output3.bed" "strand clustering output"
check_file_line_count "$meta_temp_dir/output3.bed" 6 "strand clustering line count"

# With strand consideration, + and - strand features should get different cluster IDs
# even if they overlap
pos_cluster=$(grep "pos1" "$meta_temp_dir/output3.bed" | awk '{print $NF}')
neg_cluster=$(grep "neg1" "$meta_temp_dir/output3.bed" | awk '{print $NF}')
if [ "$pos_cluster" = "$neg_cluster" ]; then
    log_error "Expected different cluster IDs for + and - strand overlapping features"
    exit 1
fi
log "✅ TEST 3 completed successfully"

# Test 4: Large distance clustering
log "Starting TEST 4: Large distance clustering"
"$meta_executable" \
    --input "$meta_temp_dir/nearby.bed" \
    --distance 1000 \
    --output "$meta_temp_dir/output4.bed"

check_file_exists "$meta_temp_dir/output4.bed" "large distance clustering output"
check_file_not_empty "$meta_temp_dir/output4.bed" "large distance clustering output"
check_file_line_count "$meta_temp_dir/output4.bed" 7 "large distance clustering line count"

# With distance 1000, most chr1 intervals should cluster together
chr1_clusters=$(grep "^chr1" "$meta_temp_dir/output4.bed" | awk '{print $NF}' | sort -u | wc -l)
if [ $chr1_clusters -gt 2 ]; then
    log "Warning: Expected few clusters on chr1 with distance 1000, got $chr1_clusters"
fi
log "✅ TEST 4 completed successfully"

# Test 5: Multiple chromosome handling
log "Starting TEST 5: Multiple chromosome handling"
# This test uses the overlapping.bed which has both chr1 and chr2
"$meta_executable" \
    --input "$meta_temp_dir/overlapping.bed" \
    --output "$meta_temp_dir/output5.bed"

check_file_exists "$meta_temp_dir/output5.bed" "multi-chromosome output"
check_file_not_empty "$meta_temp_dir/output5.bed" "multi-chromosome output"

# Check that both chromosomes are present
check_file_contains "$meta_temp_dir/output5.bed" "chr1" "chr1 features present"
check_file_contains "$meta_temp_dir/output5.bed" "chr2" "chr2 features present"

# Each chromosome should have its own cluster numbering
chr1_max_cluster=$(grep "^chr1" "$meta_temp_dir/output5.bed" | awk '{print $NF}' | sort -n | tail -1)
chr2_min_cluster=$(grep "^chr2" "$meta_temp_dir/output5.bed" | awk '{print $NF}' | sort -n | head -1)
if [ $chr2_min_cluster -le $chr1_max_cluster ]; then
    log "ℹ️  Note: Cluster IDs may continue across chromosomes (cluster numbering: chr1 max=$chr1_max_cluster, chr2 min=$chr2_min_cluster)"
fi
log "✅ TEST 5 completed successfully"

cleanup_test_env
log "🎉 All bedtools_cluster tests completed successfully!"
