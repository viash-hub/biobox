#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_flank"

# Create test data
log "Creating test data..."

# Create genome file
cat > "$meta_temp_dir/genome.txt" << 'EOF'
chr1	1000000
chr2	1000000
chr3	500000
EOF

# Create basic intervals file
cat > "$meta_temp_dir/intervals.bed" << 'EOF'
chr1	1000	2000	feature1	100	+
chr1	5000	6000	feature2	200	-
chr2	10000	11000	feature3	150	+
chr2	20000	21000	feature4	300	-
chr3	100000	101000	feature5	250	+
EOF

# Create intervals near chromosome boundaries
cat > "$meta_temp_dir/boundary.bed" << 'EOF'
chr1	10	100	start_feature	50	+
chr1	999900	999950	end_feature	75	+
chr3	490000	495000	near_end	100	+
EOF

# Create variable-sized intervals for percentage testing
cat > "$meta_temp_dir/variable.bed" << 'EOF'
chr1	10000	12000	small_2kb	10	+
chr1	20000	30000	large_10kb	20	+
chr1	50000	51000	medium_1kb	15	+
EOF

# TEST 1: Basic flanking with both sides equal
log "Starting TEST 1: Basic flanking with both sides"
"$meta_executable" \
    --input "$meta_temp_dir/intervals.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both "500" \
    --output "$meta_temp_dir/both_flanks.bed"

check_file_exists "$meta_temp_dir/both_flanks.bed" "both flanks output"
check_file_not_empty "$meta_temp_dir/both_flanks.bed" "both flanks output"

# Should create 10 intervals (5 features × 2 flanks each)
line_count=$(wc -l < "$meta_temp_dir/both_flanks.bed")
if [ "$line_count" -eq 10 ]; then
    log "✓ both flanks output has expected line count (10): $meta_temp_dir/both_flanks.bed"
else
    log "✗ both flanks output has unexpected line count ($line_count, expected 10): $meta_temp_dir/both_flanks.bed"
    exit 1
fi

log "✅ TEST 1 completed successfully"

# TEST 2: Asymmetric flanking with left and right
log "Starting TEST 2: Asymmetric flanking with left and right"
"$meta_executable" \
    --input "$meta_temp_dir/intervals.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --left "1000" \
    --right "300" \
    --output "$meta_temp_dir/asymmetric_flanks.bed"

check_file_exists "$meta_temp_dir/asymmetric_flanks.bed" "asymmetric flanks output"
check_file_not_empty "$meta_temp_dir/asymmetric_flanks.bed" "asymmetric flanks output"

# Check for different sized flanks (left flank from chr1:1000-2000 should be clamped to start at 0)
if grep -q "chr1.*0.*1000" "$meta_temp_dir/asymmetric_flanks.bed"; then
    log "✓ asymmetric flanks contains expected left flank: $meta_temp_dir/asymmetric_flanks.bed"
else
    log "✗ asymmetric flanks missing expected left flank: $meta_temp_dir/asymmetric_flanks.bed"
    cat "$meta_temp_dir/asymmetric_flanks.bed" >&2
    exit 1
fi

# Check for right flank size (300bp downstream)
if grep -q "chr1.*2000.*2300" "$meta_temp_dir/asymmetric_flanks.bed"; then
    log "✓ asymmetric flanks contains expected right flank: $meta_temp_dir/asymmetric_flanks.bed"
else
    log "✗ asymmetric flanks missing expected right flank: $meta_temp_dir/asymmetric_flanks.bed"
    cat "$meta_temp_dir/asymmetric_flanks.bed" >&2
    exit 1
fi

log "✅ TEST 2 completed successfully"

# TEST 3: Strand-aware flanking
log "Starting TEST 3: Strand-aware flanking"
"$meta_executable" \
    --input "$meta_temp_dir/intervals.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --left "800" \
    --right "400" \
    --strand \
    --output "$meta_temp_dir/strand_flanks.bed"

check_file_exists "$meta_temp_dir/strand_flanks.bed" "strand-aware flanks output"
check_file_not_empty "$meta_temp_dir/strand_flanks.bed" "strand-aware flanks output"
log "✅ TEST 3 completed successfully"

# TEST 4: Percentage-based flanking
log "Starting TEST 4: Percentage-based flanking"
"$meta_executable" \
    --input "$meta_temp_dir/variable.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both "0.5" \
    --percent \
    --output "$meta_temp_dir/percent_flanks.bed"

check_file_exists "$meta_temp_dir/percent_flanks.bed" "percentage flanks output"
check_file_not_empty "$meta_temp_dir/percent_flanks.bed" "percentage flanks output"
log "✅ TEST 4 completed successfully"

# TEST 5: Boundary handling (near chromosome ends)
log "Starting TEST 5: Boundary handling"
"$meta_executable" \
    --input "$meta_temp_dir/boundary.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both "1000" \
    --output "$meta_temp_dir/boundary_flanks.bed"

check_file_exists "$meta_temp_dir/boundary_flanks.bed" "boundary flanks output"
check_file_not_empty "$meta_temp_dir/boundary_flanks.bed" "boundary flanks output"

# Check that coordinates don't go below 0 or above chromosome length
if grep -q "^chr.*\t-" "$meta_temp_dir/boundary_flanks.bed"; then
    log "✗ boundary flanks contains negative coordinates: $meta_temp_dir/boundary_flanks.bed"
    exit 1
else
    log "✓ boundary flanks handles negative coordinates correctly: $meta_temp_dir/boundary_flanks.bed"
fi

log "✅ TEST 5 completed successfully"

# TEST 6: Header preservation
log "Starting TEST 6: Header preservation"

# Create file with header
cat > "$meta_temp_dir/with_header.bed" << 'EOF'
track name="test_track" description="Test intervals"
chr1	2000	3000	header_test	100	+
chr1	8000	9000	header_test2	150	+
EOF

"$meta_executable" \
    --input "$meta_temp_dir/with_header.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --both "200" \
    --header \
    --output "$meta_temp_dir/header_flanks.bed"

check_file_exists "$meta_temp_dir/header_flanks.bed" "header flanks output"
check_file_not_empty "$meta_temp_dir/header_flanks.bed" "header flanks output"

# Check that header is preserved
if grep -q "track name" "$meta_temp_dir/header_flanks.bed"; then
    log "✓ header flanks preserves header: $meta_temp_dir/header_flanks.bed"
else
    log "✗ header flanks missing header: $meta_temp_dir/header_flanks.bed"
    exit 1
fi

log "✅ TEST 6 completed successfully"

log "All tests completed successfully!"
