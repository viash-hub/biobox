#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_fisher"

# Create test data
log "Creating test data..."

# Create genome file
cat > "$meta_temp_dir/genome.txt" << 'EOF'
chr1	1000000
chr2	1000000
EOF

# Create file A - sorted intervals
cat > "$meta_temp_dir/intervals_a.bed" << 'EOF'
chr1	100	200	region1	10	+
chr1	300	400	region2	20	+
chr1	500	600	region3	15	-
chr2	100	200	region4	25	+
chr2	400	500	region5	30	-
EOF

# Create file B - sorted intervals with some overlaps
cat > "$meta_temp_dir/intervals_b.bed" << 'EOF'
chr1	150	250	feature1	5	+
chr1	350	450	feature2	8	+
chr1	450	550	feature3	12	-
chr2	50	150	feature4	6	+
chr2	450	550	feature5	9	-
EOF

# Create file C - larger overlap set for significance testing
cat > "$meta_temp_dir/intervals_c.bed" << 'EOF'
chr1	90	210	overlap1	10	+
chr1	290	410	overlap2	15	+
chr1	490	610	overlap3	20	-
chr2	90	210	overlap4	12	+
chr2	390	510	overlap5	18	-
chr2	600	700	overlap6	25	+
EOF

# TEST 1: Basic Fisher's exact test
log "Starting TEST 1: Basic Fisher's exact test"
"$meta_executable" \
    --input_a "$meta_temp_dir/intervals_a.bed" \
    --input_b "$meta_temp_dir/intervals_b.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/fisher_basic.txt"

check_file_exists "$meta_temp_dir/fisher_basic.txt" "basic fisher output"
check_file_not_empty "$meta_temp_dir/fisher_basic.txt" "basic fisher output"
log "✅ TEST 1 completed successfully"

# TEST 2: Fisher test with minimum overlap fraction
log "Starting TEST 2: Fisher test with overlap fractions"
"$meta_executable" \
    --input_a "$meta_temp_dir/intervals_a.bed" \
    --input_b "$meta_temp_dir/intervals_b.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --min_overlap_a 0.5 \
    --min_overlap_b 0.3 \
    --output "$meta_temp_dir/fisher_fractions.txt"

check_file_exists "$meta_temp_dir/fisher_fractions.txt" "fisher with fractions output"
check_file_not_empty "$meta_temp_dir/fisher_fractions.txt" "fisher with fractions output"
log "✅ TEST 2 completed successfully"

# TEST 3: Fisher test with reciprocal overlap
log "Starting TEST 3: Fisher test with reciprocal overlap"
"$meta_executable" \
    --input_a "$meta_temp_dir/intervals_a.bed" \
    --input_b "$meta_temp_dir/intervals_b.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --min_overlap_a 0.4 \
    --reciprocal \
    --output "$meta_temp_dir/fisher_reciprocal.txt"

check_file_exists "$meta_temp_dir/fisher_reciprocal.txt" "fisher reciprocal output"
check_file_not_empty "$meta_temp_dir/fisher_reciprocal.txt" "fisher reciprocal output"
log "✅ TEST 3 completed successfully"

# TEST 4: Fisher test with merged intervals
log "Starting TEST 4: Fisher test with merged overlapping intervals"
"$meta_executable" \
    --input_a "$meta_temp_dir/intervals_a.bed" \
    --input_b "$meta_temp_dir/intervals_c.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --merge_overlaps \
    --output "$meta_temp_dir/fisher_merged.txt"

check_file_exists "$meta_temp_dir/fisher_merged.txt" "fisher merged output"
check_file_not_empty "$meta_temp_dir/fisher_merged.txt" "fisher merged output"
log "✅ TEST 4 completed successfully"

# TEST 5: Fisher test with either overlap condition
log "Starting TEST 5: Fisher test with either overlap condition"
"$meta_executable" \
    --input_a "$meta_temp_dir/intervals_a.bed" \
    --input_b "$meta_temp_dir/intervals_b.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --min_overlap_a 0.8 \
    --min_overlap_b 0.2 \
    --either \
    --output "$meta_temp_dir/fisher_either.txt"

check_file_exists "$meta_temp_dir/fisher_either.txt" "fisher either condition output"
check_file_not_empty "$meta_temp_dir/fisher_either.txt" "fisher either condition output"
log "✅ TEST 5 completed successfully"

log "All tests completed successfully!"
