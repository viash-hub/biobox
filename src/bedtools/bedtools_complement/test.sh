#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_complement"

# Create test data
log "Creating test data..."

# Create genome file
cat > "$meta_temp_dir/genome.txt" << 'EOF'
chr1	1000
chr2	800
chr3	500
EOF

# Create simple intervals covering some regions
cat > "$meta_temp_dir/covered.bed" << 'EOF'
chr1	100	200
chr1	300	400
chr1	600	700
chr2	50	150
chr2	300	500
EOF

# Create intervals on only one chromosome
cat > "$meta_temp_dir/chr1_only.bed" << 'EOF'
chr1	100	200
chr1	500	600
chr1	800	900
EOF

# Create overlapping intervals to test merging behavior
cat > "$meta_temp_dir/overlapping.bed" << 'EOF'
chr1	100	300
chr1	250	400
chr1	600	800
chr2	100	200
chr2	150	250
EOF

# Test 1: Basic complement finding
log "Starting TEST 1: Basic complement finding"
"$meta_executable" \
    --input "$meta_temp_dir/covered.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/output1.bed"

check_file_exists "$meta_temp_dir/output1.bed" "basic complement output"
check_file_not_empty "$meta_temp_dir/output1.bed" "basic complement output"

# Should have complement regions for all chromosomes
check_file_contains "$meta_temp_dir/output1.bed" "chr1" "chr1 complement regions"
check_file_contains "$meta_temp_dir/output1.bed" "chr2" "chr2 complement regions"
check_file_contains "$meta_temp_dir/output1.bed" "chr3" "chr3 complement regions (entire chromosome)"

# Chr3 should be completely uncovered (0-500)
check_file_contains "$meta_temp_dir/output1.bed" "chr3	0	500" "complete chr3 complement"
log "✅ TEST 1 completed successfully"

# Test 2: Complement with chromosome limiting
log "Starting TEST 2: Complement with chromosome limiting"
"$meta_executable" \
    --input "$meta_temp_dir/chr1_only.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --limit_chromosomes \
    --output "$meta_temp_dir/output2.bed"

check_file_exists "$meta_temp_dir/output2.bed" "limited complement output"
check_file_not_empty "$meta_temp_dir/output2.bed" "limited complement output"

# Should only contain chr1 complement (no chr2, chr3)
check_file_contains "$meta_temp_dir/output2.bed" "chr1" "chr1 complement regions"
if grep -q "chr2\|chr3" "$meta_temp_dir/output2.bed"; then
    log_error "Expected only chr1 with -L option, but found chr2 or chr3"
    exit 1
fi
log "✅ TEST 2 completed successfully"

# Test 3: Complement of overlapping intervals
log "Starting TEST 3: Complement of overlapping intervals"
"$meta_executable" \
    --input "$meta_temp_dir/overlapping.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/output3.bed"

check_file_exists "$meta_temp_dir/output3.bed" "overlapping complement output"
check_file_not_empty "$meta_temp_dir/output3.bed" "overlapping complement output"

# bedtools complement should handle overlapping input intervals correctly
check_file_contains "$meta_temp_dir/output3.bed" "chr1" "chr1 complement with overlaps"
check_file_contains "$meta_temp_dir/output3.bed" "chr2" "chr2 complement with overlaps"
log "✅ TEST 3 completed successfully"

# Test 4: Verify complement coordinates
log "Starting TEST 4: Verify complement coordinates"
"$meta_executable" \
    --input "$meta_temp_dir/covered.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/output4.bed"

check_file_exists "$meta_temp_dir/output4.bed" "coordinate verification output"

# Check that complement starts at 0 for chr1 (nothing covered at start)
if ! grep -q "chr1	0	100" "$meta_temp_dir/output4.bed"; then
    log_error "Expected chr1 complement to start at position 0"
    exit 1
fi

# Check that complement goes to chromosome end (1000 for chr1)
if ! grep -q "700	1000" "$meta_temp_dir/output4.bed"; then
    log_error "Expected chr1 complement to end at chromosome end (1000)"
    exit 1
fi
log "✅ TEST 4 completed successfully"

# Test 5: Empty input handling
log "Starting TEST 5: Empty input handling"
# Create empty input file
touch "$meta_temp_dir/empty.bed"

"$meta_executable" \
    --input "$meta_temp_dir/empty.bed" \
    --genome "$meta_temp_dir/genome.txt" \
    --output "$meta_temp_dir/output5.bed"

check_file_exists "$meta_temp_dir/output5.bed" "empty input output"
check_file_not_empty "$meta_temp_dir/output5.bed" "empty input output"

# With no input intervals, complement should be entire genome
total_genome_size=$(awk '{sum += $2} END {print sum}' "$meta_temp_dir/genome.txt")
total_complement_size=$(awk '{sum += $3 - $2} END {print sum}' "$meta_temp_dir/output5.bed")

if [ "$total_complement_size" -ne "$total_genome_size" ]; then
    log_error "Expected complement size to equal genome size ($total_genome_size), got $total_complement_size"
    exit 1
fi
log "✅ TEST 5 completed successfully"

cleanup_test_env
log "🎉 All bedtools_complement tests completed successfully!"
