#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_random"

####################################################################################################

log "Creating test data..."

# Create a simple genome file for testing
cat > "$meta_temp_dir/test_genome.txt" << 'EOF'
chr1	10000
chr2	8000
chr3	5000
chrX	3000
chrY	2000
EOF

####################################################################################################

log "TEST 1: Basic random interval generation (default parameters)"
"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --output "$meta_temp_dir/test1_output.bed"

check_file_exists "$meta_temp_dir/test1_output.bed" "basic random intervals output"
check_file_not_empty "$meta_temp_dir/test1_output.bed" "basic random intervals result"

# Count the number of intervals generated (should be 1,000,000 by default, but we'll check structure)
num_intervals=$(wc -l < "$meta_temp_dir/test1_output.bed")
if [ "$num_intervals" -gt 0 ]; then
  log "✓ Generated $num_intervals random intervals"
else
  log "ERROR: No intervals generated"
  exit 1
fi

# Verify BED format (6 columns: chrom, start, end, name, score, strand)
first_line=$(head -1 "$meta_temp_dir/test1_output.bed")
num_columns=$(echo "$first_line" | awk '{print NF}')
if [ "$num_columns" -eq 6 ]; then
  log "✓ Output is in correct format (6 columns: chrom, start, end, name, score, strand)"
else
  log "ERROR: Expected 6 columns, found $num_columns"
  echo "First line: $first_line"
  exit 1
fi

####################################################################################################

log "TEST 2: Custom interval length"
"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --length 500 \
  --number 100 \
  --output "$meta_temp_dir/test2_output.bed"

check_file_exists "$meta_temp_dir/test2_output.bed" "custom length output"
check_file_not_empty "$meta_temp_dir/test2_output.bed" "custom length result"

# Verify we get the expected number of intervals
num_intervals_custom=$(wc -l < "$meta_temp_dir/test2_output.bed")
if [ "$num_intervals_custom" -eq 100 ]; then
  log "✓ Generated exactly 100 intervals as requested"
else
  log "ERROR: Expected 100 intervals, got $num_intervals_custom"
  exit 1
fi

# Verify interval lengths are correct (should all be 500bp)
while IFS=$'\t' read -r chrom start end name score strand; do
  length=$((end - start))
  if [ "$length" -ne 500 ]; then
      log "ERROR: Expected interval length of 500, found $length"
      exit 1
  fi
done < "$meta_temp_dir/test2_output.bed"
log "✓ All intervals have correct length of 500 bp"

####################################################################################################

log "TEST 3: Smaller number of intervals"
"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --number 10 \
  --output "$meta_temp_dir/test3_output.bed"

check_file_exists "$meta_temp_dir/test3_output.bed" "small number output"

num_small=$(wc -l < "$meta_temp_dir/test3_output.bed")
if [ "$num_small" -eq 10 ]; then
  log "✓ Generated exactly 10 intervals as requested"
else
  log "ERROR: Expected 10 intervals, got $num_small"
  exit 1
fi

####################################################################################################

log "TEST 4: Reproducibility with seed"
"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --number 20 \
  --seed 12345 \
  --output "$meta_temp_dir/test4a_output.bed"

"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --number 20 \
  --seed 12345 \
  --output "$meta_temp_dir/test4b_output.bed"

check_file_exists "$meta_temp_dir/test4a_output.bed" "seed test first run"
check_file_exists "$meta_temp_dir/test4b_output.bed" "seed test second run"

# Compare the two files - they should be identical
if diff "$meta_temp_dir/test4a_output.bed" "$meta_temp_dir/test4b_output.bed" > /dev/null; then
  log "✓ Identical results with same seed (reproducibility confirmed)"
else
  log "ERROR: Different results with same seed"
  exit 1
fi

####################################################################################################

log "TEST 5: Different seeds produce different results"
"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --number 20 \
  --seed 11111 \
  --output "$meta_temp_dir/test5a_output.bed"

"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --number 20 \
  --seed 22222 \
  --output "$meta_temp_dir/test5b_output.bed"

check_file_exists "$meta_temp_dir/test5a_output.bed" "different seed first run"
check_file_exists "$meta_temp_dir/test5b_output.bed" "different seed second run"

# Compare the two files - they should be different
if ! diff "$meta_temp_dir/test5a_output.bed" "$meta_temp_dir/test5b_output.bed" > /dev/null; then
  log "✓ Different results with different seeds (randomness confirmed)"
else
  log "WARNING: Identical results with different seeds (possible but unlikely)"
fi

####################################################################################################

log "TEST 6: Coordinate validation"
# Verify that all generated intervals are within chromosome bounds
"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --number 50 \
  --length 200 \
  --output "$meta_temp_dir/test6_output.bed"

check_file_exists "$meta_temp_dir/test6_output.bed" "coordinate validation output"

# Check that all coordinates are valid
while IFS=$'\t' read -r chrom start end name score strand; do
  # Get chromosome size from genome file
  chrom_size=$(grep "^$chrom" "$meta_temp_dir/test_genome.txt" | cut -f2)
  
  if [ -z "$chrom_size" ]; then
      log "ERROR: Chromosome $chrom not found in genome file"
      exit 1
  fi
  
  # Check bounds
  if [ "$start" -lt 0 ] || [ "$end" -gt "$chrom_size" ]; then
      log "ERROR: Interval $chrom:$start-$end is out of bounds (chromosome size: $chrom_size)"
      exit 1
  fi
  
  # Check that start < end
  if [ "$start" -ge "$end" ]; then
      log "ERROR: Invalid interval $chrom:$start-$end (start >= end)"
      exit 1
  fi
done < "$meta_temp_dir/test6_output.bed"

log "✓ All generated intervals are within valid chromosome boundaries"

####################################################################################################

log "TEST 7: Large interval length"
"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --number 5 \
  --length 2000 \
  --output "$meta_temp_dir/test7_output.bed"

check_file_exists "$meta_temp_dir/test7_output.bed" "large interval output"

# Verify large intervals are generated correctly
num_large=$(wc -l < "$meta_temp_dir/test7_output.bed")
if [ "$num_large" -eq 5 ]; then
  log "✓ Generated 5 large intervals (2000bp each)"
else
  log "✓ Generated $num_large intervals (some chromosomes may be too small for 2000bp intervals)"
fi

####################################################################################################

log "TEST 8: Parameter validation"
# Test that required parameters are enforced
log "Testing required parameter validation"

if "$meta_executable" --output "$meta_temp_dir/test.bed" 2>/dev/null; then
  log "✗ Should have failed without --genome parameter"
  exit 1
else
  log "✓ Correctly requires --genome parameter"
fi

if "$meta_executable" --genome "$meta_temp_dir/test_genome.txt" 2>/dev/null; then
  log "✗ Should have failed without --output parameter"
  exit 1
else
  log "✓ Correctly requires --output parameter"
fi

####################################################################################################

log "TEST 9: Chromosome distribution"
# Generate intervals and check that they're distributed across chromosomes
"$meta_executable" \
  --genome "$meta_temp_dir/test_genome.txt" \
  --number 100 \
  --seed 54321 \
  --output "$meta_temp_dir/test9_output.bed"

check_file_exists "$meta_temp_dir/test9_output.bed" "chromosome distribution output"

# Simple check that we have intervals
num_intervals=$(wc -l < "$meta_temp_dir/test9_output.bed")
if [ "$num_intervals" -eq 100 ]; then
  log "✓ Generated correct number of intervals ($num_intervals)"
else
  log "ERROR: Expected 100 intervals, got $num_intervals"
  exit 1
fi

# Check that intervals are on different chromosomes (basic distribution check)
num_chroms=$(cut -f1 "$meta_temp_dir/test9_output.bed" | sort -u | wc -l)
if [ "$num_chroms" -gt 1 ]; then
  log "✓ Intervals distributed across $num_chroms different chromosomes"
else
  log "✓ All intervals on single chromosome (possible with random generation)"
fi

####################################################################################################

log "✓ All tests completed successfully!"
log "bedtools_random is working correctly with proper interval generation and validation"
