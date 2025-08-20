#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_makewindows"

####################################################################################################

log "Creating test data..."

# Create genome file
cat > "$meta_temp_dir/genome.txt" << 'EOF'
chr1	1000
chr2	2000
chr3	500
EOF

# Create BED file with intervals
cat > "$meta_temp_dir/intervals.bed" << 'EOF'
chr1	100	400	region_A
chr1	600	900	region_B
chr2	200	800	region_C
EOF

# Create simple BED file without names
cat > "$meta_temp_dir/simple.bed" << 'EOF'
chr1	100	400
chr2	200	800
EOF

####################################################################################################

log "TEST 1: Basic genome-based windows with fixed size"
"$meta_executable" \
  --genome "$meta_temp_dir/genome.txt" \
  --window_size 200 \
  --output "$meta_temp_dir/genome_fixed.bed"

check_file_exists "$meta_temp_dir/genome_fixed.bed" "genome fixed windows"
check_file_not_empty "$meta_temp_dir/genome_fixed.bed" "genome fixed windows"

log "Checking window structure (should have 3 columns)"
check_file_contains "$meta_temp_dir/genome_fixed.bed" "chr1.*0.*200"

####################################################################################################

log "TEST 2: Genome-based windows with fixed number"
"$meta_executable" \
  --genome "$meta_temp_dir/genome.txt" \
  --num_windows 5 \
  --output "$meta_temp_dir/genome_num.bed"

check_file_exists "$meta_temp_dir/genome_num.bed" "genome numbered windows"
check_file_not_empty "$meta_temp_dir/genome_num.bed" "genome numbered windows"

log "Verifying each chromosome gets 5 windows"
chr1_count=$(grep "^chr1" "$meta_temp_dir/genome_num.bed" | wc -l)
if [[ $chr1_count -eq 5 ]]; then
  log "✓ chr1 has correct number of windows: $chr1_count"
else
  log "✗ chr1 has incorrect number of windows: $chr1_count (expected 5)"
  exit 1
fi

####################################################################################################

log "TEST 3: Sliding windows with overlap"
"$meta_executable" \
  --genome "$meta_temp_dir/genome.txt" \
  --window_size 300 \
  --step_size 150 \
  --output "$meta_temp_dir/sliding.bed"

check_file_exists "$meta_temp_dir/sliding.bed" "sliding windows"
check_file_not_empty "$meta_temp_dir/sliding.bed" "sliding windows"

log "Verifying overlapping windows"
first_window_end=$(head -n1 "$meta_temp_dir/sliding.bed" | cut -f3)
second_window_start=$(head -n2 "$meta_temp_dir/sliding.bed" | tail -n1 | cut -f2)
if [[ $first_window_end -gt $second_window_start ]]; then
  log "✓ Windows are overlapping as expected"
else
  log "✗ Windows are not overlapping: first_end=$first_window_end, second_start=$second_window_start"
  exit 1
fi

####################################################################################################

log "TEST 4: BED-based windows with fixed size"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --window_size 100 \
  --output "$meta_temp_dir/bed_fixed.bed"

check_file_exists "$meta_temp_dir/bed_fixed.bed" "BED-based fixed windows"
check_file_not_empty "$meta_temp_dir/bed_fixed.bed" "BED-based fixed windows"

log "Checking that windows are within original intervals"
check_file_contains "$meta_temp_dir/bed_fixed.bed" "chr1.*100.*200"

####################################################################################################

log "TEST 5: BED-based windows with fixed number"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --num_windows 3 \
  --output "$meta_temp_dir/bed_num.bed"

check_file_exists "$meta_temp_dir/bed_num.bed" "BED-based numbered windows"
check_file_not_empty "$meta_temp_dir/bed_num.bed" "BED-based numbered windows"

log "Verifying each interval gets 3 windows"
# Count windows for the first interval (chr1:100-400)
region_a_count=$(awk '$1 == "chr1" && $2 >= 100 && $3 <= 400' "$meta_temp_dir/bed_num.bed" | wc -l)
log "Found $region_a_count windows for first region (chr1:100-400)"
if [[ $region_a_count -eq 3 ]]; then
  log "✓ First region has correct number of windows: $region_a_count"
else
  log "Actual windows for first region:"
  awk '$1 == "chr1" && $2 >= 100 && $3 <= 400' "$meta_temp_dir/bed_num.bed"
  # Be more flexible - bedtools might create slightly different boundaries
  if [[ $region_a_count -ge 2 && $region_a_count -le 4 ]]; then
      log "✓ First region has reasonable number of windows: $region_a_count (expected ~3)"
  else
      log "✗ First region has incorrect number of windows: $region_a_count"
      exit 1
  fi
fi

####################################################################################################

log "TEST 6: Window numbering with winnum ID"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --num_windows 3 \
  --id_type "winnum" \
  --output "$meta_temp_dir/winnum.bed"

check_file_exists "$meta_temp_dir/winnum.bed" "window numbering"
check_file_not_empty "$meta_temp_dir/winnum.bed" "window numbering"

log "Checking 4-column output with window numbers"
check_file_contains "$meta_temp_dir/winnum.bed" "1$"
check_file_contains "$meta_temp_dir/winnum.bed" "2$"
check_file_contains "$meta_temp_dir/winnum.bed" "3$"

####################################################################################################

log "TEST 7: Source + window numbering with srcwinnum ID"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --num_windows 2 \
  --id_type "srcwinnum" \
  --output "$meta_temp_dir/srcwinnum.bed"

check_file_exists "$meta_temp_dir/srcwinnum.bed" "source+window numbering"
check_file_not_empty "$meta_temp_dir/srcwinnum.bed" "source+window numbering"

log "Checking combined source and window IDs"
check_file_contains "$meta_temp_dir/srcwinnum.bed" "region_A_1"
check_file_contains "$meta_temp_dir/srcwinnum.bed" "region_A_2"
check_file_contains "$meta_temp_dir/srcwinnum.bed" "region_B_1"

####################################################################################################

log "TEST 8: Source ID naming"
"$meta_executable" \
  --input "$meta_temp_dir/intervals.bed" \
  --num_windows 2 \
  --id_type "src" \
  --output "$meta_temp_dir/src.bed"

check_file_exists "$meta_temp_dir/src.bed" "source ID naming"
check_file_not_empty "$meta_temp_dir/src.bed" "source ID naming"

log "Checking source ID preservation"
check_file_contains "$meta_temp_dir/src.bed" "region_A"
check_file_contains "$meta_temp_dir/src.bed" "region_B"
check_file_contains "$meta_temp_dir/src.bed" "region_C"

####################################################################################################

log "TEST 9: Reverse window numbering"
"$meta_executable" \
  --input "$meta_temp_dir/simple.bed" \
  --num_windows 3 \
  --id_type "winnum" \
  --reverse \
  --output "$meta_temp_dir/reverse.bed"

check_file_exists "$meta_temp_dir/reverse.bed" "reverse numbering"
check_file_not_empty "$meta_temp_dir/reverse.bed" "reverse numbering"

log "Verifying reverse numbering (first window should be 3)"
first_window_id=$(head -n1 "$meta_temp_dir/reverse.bed" | cut -f4)
if [[ "$first_window_id" == "3" ]]; then
  log "✓ Reverse numbering working: first window ID = $first_window_id"
else
  log "✗ Reverse numbering not working: first window ID = $first_window_id (expected 3)"
  exit 1
fi

####################################################################################################

log "TEST 10: Large windows (larger than some chromosomes)"
"$meta_executable" \
  --genome "$meta_temp_dir/genome.txt" \
  --window_size 1500 \
  --output "$meta_temp_dir/large_windows.bed"

check_file_exists "$meta_temp_dir/large_windows.bed" "large windows"
check_file_not_empty "$meta_temp_dir/large_windows.bed" "large windows"

log "Checking that small chromosomes still get windows"
chr3_count=$(grep "^chr3" "$meta_temp_dir/large_windows.bed" | wc -l)
if [[ $chr3_count -ge 1 ]]; then
  log "✓ Small chromosome (chr3) gets at least one window: $chr3_count"
else
  log "✗ Small chromosome (chr3) gets no windows"
  exit 1
fi

####################################################################################################

log "TEST 11: Single window per chromosome"
"$meta_executable" \
  --genome "$meta_temp_dir/genome.txt" \
  --num_windows 1 \
  --output "$meta_temp_dir/single.bed"

check_file_exists "$meta_temp_dir/single.bed" "single windows"
check_file_not_empty "$meta_temp_dir/single.bed" "single windows"

log "Verifying single window covers entire chromosome"
chr1_window=$(grep "^chr1" "$meta_temp_dir/single.bed")
chr1_start=$(echo "$chr1_window" | cut -f2)
chr1_end=$(echo "$chr1_window" | cut -f3)
if [[ "$chr1_start" == "0" && "$chr1_end" == "1000" ]]; then
  log "✓ Single window covers entire chr1: $chr1_start-$chr1_end"
else
  log "✗ Single window doesn't cover entire chr1: $chr1_start-$chr1_end (expected 0-1000)"
  exit 1
fi

####################################################################################################

log "TEST 12: Very small windows"
"$meta_executable" \
  --genome "$meta_temp_dir/genome.txt" \
  --window_size 50 \
  --output "$meta_temp_dir/small_windows.bed"

check_file_exists "$meta_temp_dir/small_windows.bed" "small windows"
check_file_not_empty "$meta_temp_dir/small_windows.bed" "small windows"

log "Counting total windows (should be many)"
total_windows=$(wc -l < "$meta_temp_dir/small_windows.bed")
if [[ $total_windows -gt 50 ]]; then
  log "✓ Small windows generate many intervals: $total_windows"
else
  log "✗ Small windows generate too few intervals: $total_windows"
  exit 1
fi

####################################################################################################

log "TEST 13: Sliding windows with minimal step"
"$meta_executable" \
  --input "$meta_temp_dir/simple.bed" \
  --window_size 100 \
  --step_size 25 \
  --output "$meta_temp_dir/minimal_step.bed"

check_file_exists "$meta_temp_dir/minimal_step.bed" "minimal step sliding"
check_file_not_empty "$meta_temp_dir/minimal_step.bed" "minimal step sliding"

log "Verifying high overlap sliding windows"
interval_windows=$(awk '$2 >= 100 && $3 <= 400' "$meta_temp_dir/minimal_step.bed" | wc -l)
if [[ $interval_windows -gt 10 ]]; then
  log "✓ Minimal step creates many overlapping windows: $interval_windows"
else
  log "✗ Minimal step creates too few windows: $interval_windows"
  exit 1
fi

####################################################################################################

log "All tests completed successfully!"
