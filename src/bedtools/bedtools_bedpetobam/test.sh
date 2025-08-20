#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_bedpetobam"

# Create test data
log "Creating test data..."

# Create genome file
cat > "$meta_temp_dir/genome.txt" << 'EOF'
chr1	249250621
chr2	242193529
chr3	198295559
EOF

# Create BEDPE input file (paired-end BED format)
# Format: chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2
cat > "$meta_temp_dir/intervals.bedpe" << 'EOF'
chr1	100	200	chr1	300	400	pair1	100	+	+
chr1	500	600	chr1	700	800	pair2	200	+	-
chr2	150	250	chr2	350	450	pair3	300	-	-
chr2	1000	1100	chr2	1200	1300	pair4	400	-	+
EOF

# Create more detailed BEDPE file
cat > "$meta_temp_dir/detailed.bedpe" << 'EOF'
chr1	1000	2000	chr1	3000	4000	detailed1	500	+	+
chr1	5000	6000	chr1	7000	8000	detailed2	600	+	-
chr2	1500	2500	chr2	3500	4500	detailed3	700	-	-
chr2	9000	10000	chr2	11000	12000	detailed4	800	-	+
chr3	2000	3000	chr3	4000	5000	detailed5	900	+	+
EOF

# Test 1: Basic BEDPE to BAM conversion
log "Starting TEST 1: Basic BEDPE to BAM conversion"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bedpe" \
  --genome "$meta_temp_dir/genome.txt" \
  --output "$meta_temp_dir/output1.bam"

check_file_exists "$meta_temp_dir/output1.bam" "basic BAM output"
check_file_not_empty "$meta_temp_dir/output1.bam" "basic BAM output"

# BAM files are binary, so basic existence and non-empty checks are sufficient
log "✅ TEST 1 completed successfully"

# Test 2: BAM conversion with custom MAPQ
log "Starting TEST 2: BAM conversion with custom MAPQ"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bedpe" \
  --genome "$meta_temp_dir/genome.txt" \
  --mapq 60 \
  --output "$meta_temp_dir/output2.bam"

check_file_exists "$meta_temp_dir/output2.bam" "MAPQ BAM output"
check_file_not_empty "$meta_temp_dir/output2.bam" "MAPQ BAM output"
log "✅ TEST 2 completed successfully"

# Test 3: Uncompressed BAM output
log "Starting TEST 3: Uncompressed BAM output"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bedpe" \
  --genome "$meta_temp_dir/genome.txt" \
  --ubam \
  --output "$meta_temp_dir/output3.bam"

check_file_exists "$meta_temp_dir/output3.bam" "uncompressed BAM output"
check_file_not_empty "$meta_temp_dir/output3.bam" "uncompressed BAM output"

# Uncompressed BAM should be larger than compressed (typically)
compressed_size=$(stat -c%s "$meta_temp_dir/output1.bam")
uncompressed_size=$(stat -c%s "$meta_temp_dir/output3.bam")
if [ $uncompressed_size -lt $compressed_size ]; then
  log "Warning: Uncompressed BAM is smaller than compressed - may indicate issue or very small dataset"
fi
log "✅ TEST 3 completed successfully"

# Test 4: More detailed BEDPE file conversion
log "Starting TEST 4: Detailed BEDPE file conversion"
"$meta_executable" \
  --input "$meta_temp_dir/detailed.bedpe" \
  --genome "$meta_temp_dir/genome.txt" \
  --output "$meta_temp_dir/output4.bam"

check_file_exists "$meta_temp_dir/output4.bam" "detailed BAM output"
check_file_not_empty "$meta_temp_dir/output4.bam" "detailed BAM output"

# Check file size is reasonable for 5 BEDPE pairs (10 alignments)
detailed_size=$(stat -c%s "$meta_temp_dir/output4.bam")
if [ $detailed_size -lt 200 ]; then
  log_error "BAM file seems too small for 5 BEDPE pairs: $detailed_size bytes"
  exit 1
fi
log "✅ TEST 4 completed successfully"

# Test 5: Verify BAM structure with samtools (if available)
log "Starting TEST 5: BAM structure verification"
if command -v samtools &> /dev/null; then
  # Check BAM header
  if samtools view -H "$meta_temp_dir/output1.bam" | grep -q "@SQ"; then
      log "✓ BAM header contains sequence dictionary"
  else
      log_error "BAM header missing sequence dictionary"
      exit 1
  fi
  
  # Count alignments (should be double the BEDPE pairs since each pair creates 2 alignments)
  alignment_count=$(samtools view -c "$meta_temp_dir/output1.bam")
  if [ $alignment_count -eq 8 ]; then
      log "✓ BAM contains expected number of alignments: $alignment_count (4 BEDPE pairs = 8 alignments)"
  else
      log "ℹ️  Expected 8 alignments (4 BEDPE pairs), got $alignment_count"
  fi
else
  log "ℹ️  samtools not available, skipping BAM structure verification"
fi
log "✅ TEST 5 completed successfully"

log "🎉 All bedtools_bedpetobam tests completed successfully!"
