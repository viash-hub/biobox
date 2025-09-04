#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_sort"

# Create test data
log "Creating test data..."

# Create unsorted BED file for basic sorting
cat > "$meta_temp_dir/unsorted.bed" << 'EOF'
chr1	300	400
chr1	150	250  
chr1	100	200
EOF

# Create BED file with different chromosomes and sizes
cat > "$meta_temp_dir/mixed_chroms.bed" << 'EOF'
chr2	290	400
chr2	180	220
chr1	500	600
EOF

# Create BED file with scores for score-based sorting
cat > "$meta_temp_dir/with_scores.bed" << 'EOF'
chr1	100	200	feature1	960
chr1	150	250	feature2	850
chr1	300	400	feature3	740
chr2	290	390	feature4	630
chr2	180	280	feature5	920
chr3	120	220	feature6	410
EOF

# Create BED file with header
cat > "$meta_temp_dir/with_header.bed" << 'EOF'
#Header line
chr1	300	400
chr1	150	250
chr1	100	200
EOF

# Create custom genome order file
cat > "$meta_temp_dir/genome_order.txt" << 'EOF'
chr1
chr3
chr2
EOF

# Create FASTA index file
cat > "$meta_temp_dir/reference.fai" << 'EOF'
chr1	248956422
chr3	198295559
chr2	242193529
EOF

# Test 1: Basic chromosome and position sorting
log "Starting TEST 1: Basic chromosome and position sorting"
"$meta_executable" \
  --input "$meta_temp_dir/unsorted.bed" \
  --output "$meta_temp_dir/output1.bed"

check_file_exists "$meta_temp_dir/output1.bed" "basic sort output"
check_file_not_empty "$meta_temp_dir/output1.bed" "basic sort output"

# Check that features are sorted by position
check_file_contains "$meta_temp_dir/output1.bed" "chr1	100	200" "first feature by position"
check_file_line_count "$meta_temp_dir/output1.bed" 3 "basic sort line count"

# Verify the order using line numbers
head -n 1 "$meta_temp_dir/output1.bed" | grep -q "chr1	100	200" || { log "ERROR: First line is not chr1:100-200"; exit 1; }
log "✅ TEST 1 completed successfully"

# Test 2: Size-based sorting (ascending)
log "Starting TEST 2: Size-based sorting (ascending)"
"$meta_executable" \
  --input "$meta_temp_dir/mixed_chroms.bed" \
  --output "$meta_temp_dir/output2.bed" \
  --sizeA

check_file_exists "$meta_temp_dir/output2.bed" "size ascending sort output"  
check_file_not_empty "$meta_temp_dir/output2.bed" "size ascending sort output"

# Smallest feature should be first (chr2 180-220, size=40)
head -n 1 "$meta_temp_dir/output2.bed" | grep -q "chr2	180	220" || { log "ERROR: Smallest feature not first"; exit 1; }
log "✅ TEST 2 completed successfully"

# Test 3: Size-based sorting (descending)
log "Starting TEST 3: Size-based sorting (descending)"
"$meta_executable" \
  --input "$meta_temp_dir/mixed_chroms.bed" \
  --output "$meta_temp_dir/output3.bed" \
  --sizeD

check_file_exists "$meta_temp_dir/output3.bed" "size descending sort output"
check_file_not_empty "$meta_temp_dir/output3.bed" "size descending sort output"  

# Largest feature should be first (chr2 290-400, size=110)
head -n 1 "$meta_temp_dir/output3.bed" | grep -q "chr2	290	400" || { log "ERROR: Largest feature not first"; exit 1; }
log "✅ TEST 3 completed successfully"

# Test 4: Chromosome then size ascending
log "Starting TEST 4: Chromosome then size ascending"
"$meta_executable" \
  --input "$meta_temp_dir/mixed_chroms.bed" \
  --output "$meta_temp_dir/output4.bed" \
  --chrThenSizeA

check_file_exists "$meta_temp_dir/output4.bed" "chr then size asc output"
check_file_not_empty "$meta_temp_dir/output4.bed" "chr then size asc output"

# chr1 should be first, then chr2 features by size
head -n 1 "$meta_temp_dir/output4.bed" | grep -q "chr1" || { log "ERROR: chr1 not first"; exit 1; }
log "✅ TEST 4 completed successfully"

# Test 5: Score-based sorting (chromosome then score ascending)
log "Starting TEST 5: Score-based sorting (chromosome then score ascending)"
"$meta_executable" \
  --input "$meta_temp_dir/with_scores.bed" \
  --output "$meta_temp_dir/output5.bed" \
  --chrThenScoreA

check_file_exists "$meta_temp_dir/output5.bed" "chr then score asc output"
check_file_not_empty "$meta_temp_dir/output5.bed" "chr then score asc output"

# Within chr1, lowest score (740) should be first
check_file_contains "$meta_temp_dir/output5.bed" "feature3	740" "lowest score feature"
log "✅ TEST 5 completed successfully"

# Test 6: Custom genome ordering
log "Starting TEST 6: Custom genome ordering"
"$meta_executable" \
  --input "$meta_temp_dir/with_scores.bed" \
  --output "$meta_temp_dir/output6.bed" \
  --genome "$meta_temp_dir/genome_order.txt"

check_file_exists "$meta_temp_dir/output6.bed" "custom genome order output"
check_file_not_empty "$meta_temp_dir/output6.bed" "custom genome order output"

# chr1 should be first (per genome order), then chr3, then chr2
head -n 1 "$meta_temp_dir/output6.bed" | grep -q "chr1" || { log "ERROR: chr1 not first in custom order"; exit 1; }
log "✅ TEST 6 completed successfully"

# Test 7: Header preservation
log "Starting TEST 7: Header preservation"
"$meta_executable" \
  --input "$meta_temp_dir/with_header.bed" \
  --output "$meta_temp_dir/output7.bed" \
  --header

check_file_exists "$meta_temp_dir/output7.bed" "header preservation output"
check_file_not_empty "$meta_temp_dir/output7.bed" "header preservation output"

# Header should be preserved
check_file_contains "$meta_temp_dir/output7.bed" "#Header" "preserved header"
head -n 1 "$meta_temp_dir/output7.bed" | grep -q "#Header" || { log "ERROR: Header not first line"; exit 1; }
log "✅ TEST 7 completed successfully"

log "All tests completed successfully!"
