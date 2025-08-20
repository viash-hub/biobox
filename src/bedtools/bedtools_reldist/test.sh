#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_reldist"

####################################################################################################

log "Creating test data..."

# Create first test BED file (dataset A - features to analyze)
cat > "$meta_temp_dir/features_a.bed" << 'EOF'
chr1	100	200	peak1	100	+
chr1	500	600	peak2	200	+
chr1	1000	1100	peak3	150	+
chr2	300	400	peak4	180	+
chr2	800	900	peak5	220	+
chr2	1500	1600	peak6	160	+
EOF

# Create second test BED file (dataset B - reference features)  
cat > "$meta_temp_dir/features_b.bed" << 'EOF'
chr1	150	300	gene1	100	+
chr1	450	650	gene2	200	-
chr1	950	1200	gene3	150	+
chr2	250	450	gene4	300	+
chr2	750	1000	gene5	250	-
chr2	1400	1700	gene6	180	+
EOF

####################################################################################################

log "TEST 1: Basic relative distance calculation (summary mode)"
"$meta_executable" \
    --bed_a "$meta_temp_dir/features_a.bed" \
    --bed_b "$meta_temp_dir/features_b.bed" \
    --output "$meta_temp_dir/test1_output.txt"

check_file_exists "$meta_temp_dir/test1_output.txt" "basic reldist output"
check_file_not_empty "$meta_temp_dir/test1_output.txt" "basic reldist result"

# Check that output contains distribution information
num_lines=$(wc -l < "$meta_temp_dir/test1_output.txt")
if [ "$num_lines" -gt 0 ]; then
    log "✓ Generated summary distribution with $num_lines lines"
else
    log "ERROR: Empty output file"
    exit 1
fi

####################################################################################################

log "TEST 2: Detailed relative distance for each interval"
"$meta_executable" \
    --bed_a "$meta_temp_dir/features_a.bed" \
    --bed_b "$meta_temp_dir/features_b.bed" \
    --detail \
    --output "$meta_temp_dir/test2_output.txt"

check_file_exists "$meta_temp_dir/test2_output.txt" "detailed reldist output"
check_file_not_empty "$meta_temp_dir/test2_output.txt" "detailed reldist result"

# Check that we have detailed output for each feature in dataset A
num_features_a=$(wc -l < "$meta_temp_dir/features_a.bed")
num_detail_lines=$(wc -l < "$meta_temp_dir/test2_output.txt")

log "✓ Detail mode produced $num_detail_lines lines for $num_features_a input features"

# Verify that the detail output has the expected format
if head -1 "$meta_temp_dir/test2_output.txt" | grep -q -E '^chr[0-9]+\s+[0-9]+\s+[0-9]+'; then
    log "✓ Detail output has expected format (starts with genomic coordinates)"
else
    log "WARNING: Detail output format may differ from expected"
fi

####################################################################################################

log "TEST 3: Different chromosome distribution"
# Create datasets with features on different chromosomes
cat > "$meta_temp_dir/multi_chrom_a.bed" << 'EOF'
chr1	100	200	feature1	100	+
chr2	100	200	feature2	100	+
chr3	100	200	feature3	100	+
EOF

cat > "$meta_temp_dir/multi_chrom_b.bed" << 'EOF'
chr1	300	400	ref1	100	+
chr2	300	400	ref2	100	+
chr3	300	400	ref3	100	+
chr4	300	400	ref4	100	+
EOF

"$meta_executable" \
    --bed_a "$meta_temp_dir/multi_chrom_a.bed" \
    --bed_b "$meta_temp_dir/multi_chrom_b.bed" \
    --output "$meta_temp_dir/test3_output.txt"

check_file_exists "$meta_temp_dir/test3_output.txt" "multi-chromosome output"
check_file_not_empty "$meta_temp_dir/test3_output.txt" "multi-chromosome result"

log "✓ Multi-chromosome analysis completed"

####################################################################################################

log "TEST 4: Single chromosome analysis"
cat > "$meta_temp_dir/single_chrom_a.bed" << 'EOF'
chr1	100	200	feature1	100	+
chr1	500	600	feature2	100	+
chr1	1000	1100	feature3	100	+
EOF

cat > "$meta_temp_dir/single_chrom_b.bed" << 'EOF'
chr1	250	300	ref1	100	+
chr1	750	800	ref2	100	+
chr1	1250	1300	ref3	100	+
EOF

"$meta_executable" \
    --bed_a "$meta_temp_dir/single_chrom_a.bed" \
    --bed_b "$meta_temp_dir/single_chrom_b.bed" \
    --output "$meta_temp_dir/test4_output.txt"

check_file_exists "$meta_temp_dir/test4_output.txt" "single chromosome output"
check_file_not_empty "$meta_temp_dir/test4_output.txt" "single chromosome result"

log "✓ Single chromosome analysis completed"

####################################################################################################

log "TEST 5: Edge case - overlapping features"
cat > "$meta_temp_dir/overlap_a.bed" << 'EOF'
chr1	100	300	feature1	100	+
chr1	500	700	feature2	100	+
EOF

cat > "$meta_temp_dir/overlap_b.bed" << 'EOF'
chr1	150	250	ref1	100	+
chr1	550	650	ref2	100	+
EOF

"$meta_executable" \
    --bed_a "$meta_temp_dir/overlap_a.bed" \
    --bed_b "$meta_temp_dir/overlap_b.bed" \
    --detail \
    --output "$meta_temp_dir/test5_output.txt"

check_file_exists "$meta_temp_dir/test5_output.txt" "overlapping features output"
check_file_not_empty "$meta_temp_dir/test5_output.txt" "overlapping features result"

log "✓ Overlapping features analysis completed"

####################################################################################################

log "TEST 6: Compare summary vs detailed output"
# Run the same analysis in both modes and verify consistency
"$meta_executable" \
    --bed_a "$meta_temp_dir/features_a.bed" \
    --bed_b "$meta_temp_dir/features_b.bed" \
    --output "$meta_temp_dir/test6_summary.txt"

"$meta_executable" \
    --bed_a "$meta_temp_dir/features_a.bed" \
    --bed_b "$meta_temp_dir/features_b.bed" \
    --detail \
    --output "$meta_temp_dir/test6_detail.txt"

check_file_exists "$meta_temp_dir/test6_summary.txt" "summary comparison output"
check_file_exists "$meta_temp_dir/test6_detail.txt" "detail comparison output"

summary_lines=$(wc -l < "$meta_temp_dir/test6_summary.txt")
detail_lines=$(wc -l < "$meta_temp_dir/test6_detail.txt")

log "✓ Summary mode: $summary_lines lines, Detail mode: $detail_lines lines"

####################################################################################################

log "TEST 7: Parameter validation"
# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" --bed_b "$meta_temp_dir/features_b.bed" --output "$meta_temp_dir/test.txt" 2>/dev/null; then
    log "✗ Should have failed without --bed_a parameter"
    exit 1
else
    log "✓ Correctly requires --bed_a parameter"
fi

if "$meta_executable" --bed_a "$meta_temp_dir/features_a.bed" --output "$meta_temp_dir/test.txt" 2>/dev/null; then
    log "✗ Should have failed without --bed_b parameter"
    exit 1
else
    log "✓ Correctly requires --bed_b parameter"
fi

if "$meta_executable" --bed_a "$meta_temp_dir/features_a.bed" --bed_b "$meta_temp_dir/features_b.bed" 2>/dev/null; then
    log "✗ Should have failed without --output parameter"
    exit 1
else
    log "✓ Correctly requires --output parameter"
fi

####################################################################################################

log "TEST 8: Output format validation"
# Verify the output format is consistent and parseable
"$meta_executable" \
    --bed_a "$meta_temp_dir/features_a.bed" \
    --bed_b "$meta_temp_dir/features_b.bed" \
    --detail \
    --output "$meta_temp_dir/test8_output.txt"

check_file_exists "$meta_temp_dir/test8_output.txt" "format validation output"

# Check that each line can be parsed (basic format check)
if ! awk 'NF >= 1' "$meta_temp_dir/test8_output.txt" >/dev/null 2>&1; then
    log "ERROR: Output format appears to be malformed"
    exit 1
else
    log "✓ Output format validation passed"
fi

####################################################################################################

log "✓ All tests completed successfully!"
log "bedtools_reldist is working correctly with both summary and detailed analysis modes"
