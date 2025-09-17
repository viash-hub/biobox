#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_merge"

# Create test data
log "Creating test data..."

# Create basic BED file with overlapping features
cat > "$meta_temp_dir/featureA.bed" << 'EOF'
chr1	100	200
chr1	150	250
chr1	300	400
EOF

# Create BED file with strand information  
cat > "$meta_temp_dir/featureB.bed" << 'EOF'
chr1	100	200	a1	1	+
chr1	180	250	a2	2	+
chr1	250	500	a3	3	-
chr1	501	1000	a4	4	+
EOF

# Create BED file for precision testing
cat > "$meta_temp_dir/feature_precision.bed" << 'EOF'
chr1	100	200	a1	1.9	+
chr1	180	250	a2	2.5	+
chr1	250	500	a3	3.3	-
chr1	501	1000	a4	4	+
EOF

# Create GFF file for header testing
cat > "$meta_temp_dir/feature.gff" << 'EOF'
##gff-version 3
chr1	.	gene	1000	2000	.	+	.	ID=gene1;Name=Gene1
chr1	.	exon	1000	1200	.	+	.	ID=exon1;Parent=transcript1
chr1	.	CDS	1000	1200	.	+	0	ID=cds1;Parent=transcript1
chr1	.	CDS	1500	1700	.	+	2	ID=cds2;Parent=transcript1
chr2	.	exon	1500	1700	.	+	.	ID=exon2;Parent=transcript1
chr3	.	mRNA	1000	2000	.	+	.	ID=transcript1;Parent=gene1
EOF

# Test 1: Basic merge functionality
log "Starting TEST 1: Basic merge functionality"
"$meta_executable" \
  --input "$meta_temp_dir/featureA.bed" \
  --output "$meta_temp_dir/output1.bed"

check_file_exists "$meta_temp_dir/output1.bed" "basic merge output"
check_file_not_empty "$meta_temp_dir/output1.bed" "basic merge output"

# Check that the merged intervals are as expected
check_file_contains "$meta_temp_dir/output1.bed" "chr1	100	250" "first merged interval"
check_file_contains "$meta_temp_dir/output1.bed" "chr1	300	400" "second merged interval"
check_file_line_count "$meta_temp_dir/output1.bed" 2 "merged output line count"
log "✅ TEST 1 completed successfully"

# Test 2: Strand-specific merging
log "Starting TEST 2: Strand-specific merging"
"$meta_executable" \
  --input "$meta_temp_dir/featureB.bed" \
  --output "$meta_temp_dir/output2.bed" \
  --strand

check_file_exists "$meta_temp_dir/output2.bed" "strand-specific merge output"
check_file_not_empty "$meta_temp_dir/output2.bed" "strand-specific merge output"

# Check that strand-specific merging occurred
check_file_contains "$meta_temp_dir/output2.bed" "chr1	100	250" "merged + strand features"
check_file_contains "$meta_temp_dir/output2.bed" "chr1	250	500" "- strand feature"
check_file_contains "$meta_temp_dir/output2.bed" "chr1	501	1000" "+ strand feature"
log "✅ TEST 2 completed successfully"

# Test 3: Distance-based merging
log "Starting TEST 3: Distance-based merging"  
"$meta_executable" \
  --input "$meta_temp_dir/featureA.bed" \
  --output "$meta_temp_dir/output3.bed" \
  --distance 50

check_file_exists "$meta_temp_dir/output3.bed" "distance-based merge output"
check_file_not_empty "$meta_temp_dir/output3.bed" "distance-based merge output"

# Expected: all features merged into one (distance allows gap between 250-300)
check_file_contains "$meta_temp_dir/output3.bed" "chr1	100	400" "distance-based merged interval"
check_file_line_count "$meta_temp_dir/output3.bed" 1 "distance-based merge line count"
log "✅ TEST 3 completed successfully"

# Test 4: Column operations with aggregation
log "Starting TEST 4: Column operations with aggregation"
"$meta_executable" \
  --input "$meta_temp_dir/featureB.bed" \
  --output "$meta_temp_dir/output4.bed" \
  --columns "5" \
  --operation "mean"

check_file_exists "$meta_temp_dir/output4.bed" "column aggregation output"
check_file_not_empty "$meta_temp_dir/output4.bed" "column aggregation output"

# Check that output contains numerical values (mean of column 5)
check_file_contains "$meta_temp_dir/output4.bed" "1.5" "mean aggregation result"
log "✅ TEST 4 completed successfully"

# Test 5: Custom delimiter for collapse operation
log "Starting TEST 5: Custom delimiter for collapse operation"
"$meta_executable" \
  --input "$meta_temp_dir/featureB.bed" \
  --output "$meta_temp_dir/output5.bed" \
  --columns "4" \
  --operation "collapse" \
  --delimiter "|"

check_file_exists "$meta_temp_dir/output5.bed" "custom delimiter output"
check_file_not_empty "$meta_temp_dir/output5.bed" "custom delimiter output"

# Check that output contains pipe-separated values
check_file_contains "$meta_temp_dir/output5.bed" "|" "pipe delimiter in collapsed values"
log "✅ TEST 5 completed successfully"

log "All tests completed successfully!"
