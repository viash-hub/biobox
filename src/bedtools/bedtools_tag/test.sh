#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_tag"

####################################################################################################

log "Creating test data..."

# Create a simple BED file representing aligned reads
log "Creating test aligned reads as BED..."
cat > "$meta_temp_dir/alignments.bed" << 'EOF'
chr22	1100	1200	read1	60	+
chr22	1800	1900	read2	60	-
chr22	3200	3300	read3	60	+
chr22	4100	4200	read4	60	-
chr22	5200	5300	read5	60	+
chr22	6100	6200	read6	60	-
EOF

# Create a genome file for bedToBam
cat > "$meta_temp_dir/genome.txt" << 'EOF'
chr22	50000000
EOF

# Convert BED to BAM using bedtools
log "Converting BED to BAM..."
bedToBam -i "$meta_temp_dir/alignments.bed" -g "$meta_temp_dir/genome.txt" > "$meta_temp_dir/input.bam"

# Create annotation BED files that overlap with our test reads
cat > "$meta_temp_dir/genes.bed" << 'EOF'
chr22	1050	1250	gene1	100	+
chr22	3150	3350	gene2	200	+
chr22	5150	5350	gene3	300	-
EOF

cat > "$meta_temp_dir/enhancers.bed" << 'EOF'
chr22	1750	1950	enhancer1	150	+
chr22	4050	4250	enhancer2	250	+
chr22	6050	6250	enhancer3	350	+
EOF

check_file_exists "$meta_temp_dir/input.bam" "input BAM file"

####################################################################################################

log "TEST 1: Basic tagging with single BED file"

"$meta_executable" \
  --input "$meta_temp_dir/input.bam" \
  --files "$meta_temp_dir/genes.bed" \
  --labels "GENE" \
  --output "$meta_temp_dir/tagged_basic.bam"

check_file_exists "$meta_temp_dir/tagged_basic.bam" "basic tagged output"
check_file_not_empty "$meta_temp_dir/tagged_basic.bam" "basic tagged output"

log "✓ Basic tagging test passed"

####################################################################################################

log "TEST 2: Custom tag name"

"$meta_executable" \
  --input "$meta_temp_dir/input.bam" \
  --files "$meta_temp_dir/genes.bed" \
  --labels "GENE" \
  --output "$meta_temp_dir/tagged_custom.bam" \
  --tag_name "XG"

check_file_exists "$meta_temp_dir/tagged_custom.bam" "custom tag output"
check_file_not_empty "$meta_temp_dir/tagged_custom.bam" "custom tag output"

log "✓ Custom tag name test passed"

####################################################################################################

log "TEST 3: Multiple annotation files with names"

"$meta_executable" \
  --input "$meta_temp_dir/input.bam" \
  --files "$meta_temp_dir/genes.bed;$meta_temp_dir/enhancers.bed" \
  --output "$meta_temp_dir/tagged_names.bam" \
  --use_names

check_file_exists "$meta_temp_dir/tagged_names.bam" "names tagged output"
check_file_not_empty "$meta_temp_dir/tagged_names.bam" "names tagged output"

log "✓ Multiple files with names test passed"

####################################################################################################

log "TEST 4: Include scores option"

"$meta_executable" \
  --input "$meta_temp_dir/input.bam" \
  --files "$meta_temp_dir/genes.bed" \
  --output "$meta_temp_dir/tagged_scores.bam" \
  --use_scores

check_file_exists "$meta_temp_dir/tagged_scores.bam" "scores tagged output"
check_file_not_empty "$meta_temp_dir/tagged_scores.bam" "scores tagged output"

log "✓ Include scores test passed"

####################################################################################################

log "TEST 5: Overlap mode testing"

"$meta_executable" \
  --input "$meta_temp_dir/input.bam" \
  --files "$meta_temp_dir/genes.bed" \
  --labels "GENE" \
  --output "$meta_temp_dir/tagged_overlap.bam" \
  --min_overlap 0.5

check_file_exists "$meta_temp_dir/tagged_overlap.bam" "overlap tagged output"
check_file_not_empty "$meta_temp_dir/tagged_overlap.bam" "overlap tagged output"

log "✓ Overlap mode test passed"

####################################################################################################

log "TEST 6: Error handling - Missing input file"

if "$meta_executable" \
  --input "/nonexistent/file.bam" \
  --files "$meta_temp_dir/genes.bed" \
  --labels "GENE" \
  --output "$meta_temp_dir/error_test.bam" 2>/dev/null; then
  log "✗ Should have failed with missing input file"
  exit 1
else
  log "✓ Correctly handled missing input file"
fi

####################################################################################################

log "TEST 7: Error handling - Missing annotation file"

if "$meta_executable" \
  --input "$meta_temp_dir/input.bam" \
  --files "/nonexistent/file.bed" \
  --labels "GENE" \
  --output "$meta_temp_dir/error_test.bam" 2>/dev/null; then
  log "✗ Should have failed with missing annotation file"
  exit 1
else
  log "✓ Correctly handled missing annotation file"
fi

####################################################################################################

log "TEST 8: Error handling - Missing output parameter"

if "$meta_executable" \
  --input "$meta_temp_dir/input.bam" \
  --files "$meta_temp_dir/genes.bed" \
  --labels "GENE" 2>/dev/null; then
  log "✗ Should have failed without output parameter"
  exit 1
else
  log "✓ Correctly handled missing output parameter"
fi

####################################################################################################

log "All tests completed successfully!"
