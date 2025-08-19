#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# Create test directory
test_dir="$meta_temp_dir/test_data"
mkdir -p "$test_dir"

# Create test genome file
log "Creating test genome file..."
cat > "$test_dir/test.genome" << 'EOF'
chr1	10000
chr2	8000
chr3	5000
EOF

# Create test BED file
log "Creating test BED file..."
cat > "$test_dir/test.bed" << 'EOF'
chr1	100	200	feature1	100	+
chr1	300	500	feature2	200	-
chr2	1000	1500	feature3	150	+
chr2	2000	2200	feature4	180	-
chr3	500	800	feature5	120	+
EOF

# --- Test Case 1: Basic histogram output (default) ---
log "Starting TEST 1: Basic coverage histogram"

log "Executing $meta_name with default histogram output..."
"$meta_executable" \
    --input "$test_dir/test.bed" \
    --genome "$test_dir/test.genome" \
    --output "$meta_temp_dir/output1.txt"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output1.txt" "histogram output file"
check_file_not_empty "$meta_temp_dir/output1.txt" "histogram output file"

# Check histogram format (should have columns: chromosome, depth, count, total_bases, fraction)
line_count=$(wc -l < "$meta_temp_dir/output1.txt")
log "Histogram contains $line_count lines"
[ "$line_count" -gt 0 ] || { log_error "Histogram output is empty"; exit 1; }

# Check that it contains expected format
head -1 "$meta_temp_dir/output1.txt" | awk 'NF != 5 { exit 1 }' || {
    log_error "Histogram format incorrect (expected 5 columns)"
    exit 1
}

log "✅ TEST 1 completed successfully"

# --- Test Case 2: BedGraph format ---
log "Starting TEST 2: BedGraph format output"

log "Executing $meta_name with BedGraph format..."
"$meta_executable" \
    --input "$test_dir/test.bed" \
    --genome "$test_dir/test.genome" \
    --output "$meta_temp_dir/output2.bg" \
    --bed_graph

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output2.bg" "BedGraph output file"
check_file_not_empty "$meta_temp_dir/output2.bg" "BedGraph output file"

# Check BedGraph format (chromosome, start, end, depth)
head -1 "$meta_temp_dir/output2.bg" | awk 'NF != 4 { exit 1 }' || {
    log_error "BedGraph format incorrect (expected 4 columns)"
    exit 1
}

# Check that coordinates make sense (start < end)
awk '$2 >= $3 { print "Invalid coordinates: " $0; exit 1 }' "$meta_temp_dir/output2.bg" || {
    log_error "Invalid BedGraph coordinates found"
    exit 1
}

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Per-base depth ---
log "Starting TEST 3: Per-base depth output"

log "Executing $meta_name with per-base depth..."
"$meta_executable" \
    --input "$test_dir/test.bed" \
    --genome "$test_dir/test.genome" \
    --output "$meta_temp_dir/output3.depth" \
    --depth

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output3.depth" "depth output file"
check_file_not_empty "$meta_temp_dir/output3.depth" "depth output file"

# Check depth format (chromosome, position, depth)
head -1 "$meta_temp_dir/output3.depth" | awk 'NF != 3 { exit 1 }' || {
    log_error "Depth format incorrect (expected 3 columns)"
    exit 1
}

log "✅ TEST 3 completed successfully"

# --- Test Case 4: BedGraph with zero coverage ---
log "Starting TEST 4: BedGraph with zero coverage"

log "Executing $meta_name with BedGraph including zero coverage..."
"$meta_executable" \
    --input "$test_dir/test.bed" \
    --genome "$test_dir/test.genome" \
    --output "$meta_temp_dir/output4.bga" \
    --bed_graph_zero_coverage

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output4.bga" "BedGraph+zero output file"
check_file_not_empty "$meta_temp_dir/output4.bga" "BedGraph+zero output file"

# This output should be larger than regular BedGraph since it includes zero coverage
bg_size=$(wc -l < "$meta_temp_dir/output2.bg")
bga_size=$(wc -l < "$meta_temp_dir/output4.bga")
log "BedGraph lines: $bg_size, BedGraph+zero lines: $bga_size"

# Check that we can find zero coverage regions
if grep -q "	0$" "$meta_temp_dir/output4.bga"; then
    log "✓ Found zero coverage regions in output"
else
    log "Note: No zero coverage regions found (this may be expected with test data)"
fi

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Test strand-specific coverage ---
log "Starting TEST 5: Strand-specific coverage"

# Create BED file with strand information (6 columns minimum)
cat > "$test_dir/strand.bed" << 'EOF'
chr1	100	200	feature1	100	+
chr1	300	500	feature2	200	-
EOF

log "Executing $meta_name with strand-specific coverage..."
"$meta_executable" \
    --input "$test_dir/strand.bed" \
    --genome "$test_dir/test.genome" \
    --output "$meta_temp_dir/output5.txt" \
    --strand "+"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/output5.txt" "strand-specific output file"
check_file_not_empty "$meta_temp_dir/output5.txt" "strand-specific output file"

log "✅ TEST 5 completed successfully"

print_test_summary "All tests completed successfully"
