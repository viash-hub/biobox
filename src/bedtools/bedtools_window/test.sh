#!/bin/bash

## VIASH START
## VIASH END

# Load test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
init_test_env

log "Starting tests for bedtools_window"

# Create test data
log "Creating test data..."

# Create query file A with features to create windows around
cat > "$meta_temp_dir/queries.bed" << 'EOF'
chr1	200	300	query1	100	+
chr1	500	600	query2	200	+
chr1	800	900	query3	300	-
chr2	100	200	query4	150	-
EOF

# Create database file B with features to find within windows
cat > "$meta_temp_dir/features.bed" << 'EOF'
chr1	150	250	feature1	50	+
chr1	180	220	feature2	60	-
chr1	400	450	feature3	70	+
chr1	550	650	feature4	80	-
chr1	750	850	feature5	90	+
chr1	820	920	feature6	100	-
chr2	50	120	feature7	110	+
chr2	150	250	feature8	120	-
EOF

####################################################################################################

log "TEST 1: Basic window search with default settings"

"$meta_executable" \
  --input_a "$meta_temp_dir/queries.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --output "$meta_temp_dir/output1.bed"

check_file_exists "$meta_temp_dir/output1.bed" "basic window output"
check_file_not_empty "$meta_temp_dir/output1.bed" "basic window output"

# Should find overlaps within default 1000bp windows
line_count=$(wc -l < "$meta_temp_dir/output1.bed")
if [ "$line_count" -lt 1 ]; then
  log "❌ Expected at least 1 overlap with default 1000bp windows"
  exit 1
fi

log "✅ TEST 1 completed successfully"

####################################################################################################

log "TEST 2: Custom symmetric window size"

"$meta_executable" \
  --input_a "$meta_temp_dir/queries.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --window_size 100 \
  --output "$meta_temp_dir/output2.bed"

check_file_exists "$meta_temp_dir/output2.bed" "custom window output"
check_file_not_empty "$meta_temp_dir/output2.bed" "custom window output"

# With 100bp windows, should have fewer overlaps than default
line_count_small=$(wc -l < "$meta_temp_dir/output2.bed")
line_count_default=$(wc -l < "$meta_temp_dir/output1.bed")

log "✅ TEST 2 completed successfully"

####################################################################################################

log "TEST 3: Asymmetric windows"

"$meta_executable" \
  --input_a "$meta_temp_dir/queries.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --left_window 50 \
  --right_window 200 \
  --output "$meta_temp_dir/output3.bed"

check_file_exists "$meta_temp_dir/output3.bed" "asymmetric window output"
check_file_not_empty "$meta_temp_dir/output3.bed" "asymmetric window output"

log "✅ TEST 3 completed successfully"

####################################################################################################

log "TEST 4: Strand-specific overlaps (same strand)"

"$meta_executable" \
  --input_a "$meta_temp_dir/queries.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --same_strand \
  --output "$meta_temp_dir/output4.bed"

check_file_exists "$meta_temp_dir/output4.bed" "same strand output"

# Check that all overlaps have matching strand (6th and 12th columns should match)
if [ -s "$meta_temp_dir/output4.bed" ]; then
  while IFS=$'\t' read -r c1 s1 e1 n1 sc1 str1 c2 s2 e2 n2 sc2 str2 rest; do
      if [ "$str1" != "$str2" ]; then
      log "❌ Found opposite strand overlap in same-strand mode: $str1 vs $str2"
      exit 1
      fi
  done < "$meta_temp_dir/output4.bed"
fi

log "✅ TEST 4 completed successfully"

####################################################################################################

log "TEST 5: Count overlaps mode"

"$meta_executable" \
  --input_a "$meta_temp_dir/queries.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --count \
  --output "$meta_temp_dir/output5.bed"

check_file_exists "$meta_temp_dir/output5.bed" "count mode output"
check_file_not_empty "$meta_temp_dir/output5.bed" "count mode output"

# Should have exactly 4 lines (one per query)
check_file_line_count "$meta_temp_dir/output5.bed" 4 "count mode line count"

# Check that each line has a count column (7th column should be numeric)
while IFS=$'\t' read -r c s e n sc str count rest; do
  if ! [[ "$count" =~ ^[0-9]+$ ]]; then
      log "❌ Count column should be numeric, got: $count"
      exit 1
  fi
done < "$meta_temp_dir/output5.bed"

log "✅ TEST 5 completed successfully"

####################################################################################################

log "TEST 6: Unique mode (report A entries once)"

"$meta_executable" \
  --input_a "$meta_temp_dir/queries.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --unique \
  --output "$meta_temp_dir/output6.bed"

check_file_exists "$meta_temp_dir/output6.bed" "unique mode output"

# Should have at most 4 lines (one per query that has overlaps)
line_count=$(wc -l < "$meta_temp_dir/output6.bed")
if [ "$line_count" -gt 4 ]; then
  log "❌ Unique mode should have at most 4 lines, got $line_count"
  exit 1
fi

# Each line should have only A file columns (6 columns)
if [ "$line_count" -gt 0 ]; then
  col_count=$(head -1 "$meta_temp_dir/output6.bed" | awk -F'\t' '{print NF}')
  if [ "$col_count" -ne 6 ]; then
      log "❌ Unique mode should have 6 columns (A file only), got $col_count"
      exit 1
  fi
fi

log "✅ TEST 6 completed successfully"

####################################################################################################

log "TEST 7: No overlaps mode (features with no nearby features)"

"$meta_executable" \
  --input_a "$meta_temp_dir/queries.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --no_overlaps \
  --window_size 10 \
  --output "$meta_temp_dir/output7.bed"

check_file_exists "$meta_temp_dir/output7.bed" "no overlaps output"

# With very small windows (10bp), there should be some queries with no overlaps
# Each line should have only A file columns (6 columns)
if [ -s "$meta_temp_dir/output7.bed" ]; then
  col_count=$(head -1 "$meta_temp_dir/output7.bed" | awk -F'\t' '{print NF}')
  if [ "$col_count" -ne 6 ]; then
      log "❌ No overlaps mode should have 6 columns (A file only), got $col_count"
      exit 1
  fi
fi

log "✅ TEST 7 completed successfully"

####################################################################################################

log "TEST 8: Header preservation"

# Create input file with header
cat > "$meta_temp_dir/queries_with_header.bed" << 'EOF'
#chrom	start	end	name	score	strand
chr1	200	300	query1	100	+
chr1	500	600	query2	200	+
EOF

"$meta_executable" \
  --input_a "$meta_temp_dir/queries_with_header.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --header \
  --output "$meta_temp_dir/output8.bed"

check_file_exists "$meta_temp_dir/output8.bed" "header output"
check_file_not_empty "$meta_temp_dir/output8.bed" "header output"

# First line should be header (start with #)
first_line=$(head -1 "$meta_temp_dir/output8.bed")
if [[ ! "$first_line" =~ ^#.* ]]; then
  log "❌ Header should start with #, got: $first_line"
  exit 1
fi

log "✅ TEST 8 completed successfully"

####################################################################################################

log "TEST 9: Error handling - Missing input file"

if "$meta_executable" \
  --input_a "/nonexistent/file.bed" \
  --input_b "$meta_temp_dir/features.bed" \
  --output "$meta_temp_dir/error_test.bed" 2>/dev/null; then
  log "❌ Should fail with missing input file"
  exit 1
fi

log "✓ Correctly handled missing input file"

####################################################################################################

log "TEST 10: Error handling - Missing output parameter"

if "$meta_executable" \
  --input_a "$meta_temp_dir/queries.bed" \
  --input_b "$meta_temp_dir/features.bed" 2>/dev/null; then
  log "❌ Should fail without output parameter"
  exit 1
fi

log "✓ Correctly handled missing output parameter"

log "All tests completed successfully!"
