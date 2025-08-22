#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

log "Starting tests for $meta_name"

####################################################################################################

log "Creating test data..."

# Create test BEDGRAPH files with overlapping and non-overlapping intervals
cat > "$meta_temp_dir/sample1.bedgraph" << 'EOF'
chr1	100	200	1.5
chr1	300	400	2.0
chr1	500	600	3.5
EOF

cat > "$meta_temp_dir/sample2.bedgraph" << 'EOF'
chr1	150	250	4.0
chr1	350	450	5.5
chr1	550	650	6.0
EOF

cat > "$meta_temp_dir/sample3.bedgraph" << 'EOF'
chr1	175	275	7.0
chr1	375	475	8.5
chr1	575	675	9.0
EOF

####################################################################################################

log "TEST 1: Basic union of two BEDGRAPH files"

"$meta_executable" \
  --files "$meta_temp_dir/sample1.bedgraph;$meta_temp_dir/sample2.bedgraph" \
  --output "$meta_temp_dir/union_basic.bedgraph"

check_file_exists "$meta_temp_dir/union_basic.bedgraph" "basic union output"
check_file_not_empty "$meta_temp_dir/union_basic.bedgraph" "basic union output"

log "✓ Basic union test passed"

####################################################################################################

log "TEST 2: Union with header option"

"$meta_executable" \
  --files "$meta_temp_dir/sample1.bedgraph;$meta_temp_dir/sample2.bedgraph" \
  --output "$meta_temp_dir/union_header.bedgraph" \
  --header

check_file_exists "$meta_temp_dir/union_header.bedgraph" "union output with header"
check_file_not_empty "$meta_temp_dir/union_header.bedgraph" "union output with header"

# Check that header line is present (should start with 'chrom')
if head -1 "$meta_temp_dir/union_header.bedgraph" | grep -q "chrom"; then
  log "✓ Header line correctly added"
else
  log "✗ Header line missing or incorrect"
  exit 1
fi

log "✓ Header option test passed"

####################################################################################################

log "TEST 3: Union with multiple files"

"$meta_executable" \
  --files "$meta_temp_dir/sample1.bedgraph;$meta_temp_dir/sample2.bedgraph;$meta_temp_dir/sample3.bedgraph" \
  --output "$meta_temp_dir/union_multiple.bedgraph"

check_file_exists "$meta_temp_dir/union_multiple.bedgraph" "multiple files union output"
check_file_not_empty "$meta_temp_dir/union_multiple.bedgraph" "multiple files union output"

# Check that output has the expected number of columns (3 for coordinates + 3 for values)
expected_columns=6
actual_columns=$(head -1 "$meta_temp_dir/union_multiple.bedgraph" | wc -w)
if [ "$actual_columns" -eq "$expected_columns" ]; then
  log "✓ Output has correct number of columns ($actual_columns)"
else
  log "✗ Expected $expected_columns columns, got $actual_columns"
  exit 1
fi

log "✓ Multiple files test passed"

####################################################################################################

log "TEST 4: Error handling - Missing input files"

if "$meta_executable" \
  --files "/nonexistent/file1.bedgraph" "/nonexistent/file2.bedgraph" \
  --output "$meta_temp_dir/error_test.bedgraph" 2>/dev/null; then
  log "✗ Should have failed with missing input files"
  exit 1
else
  log "✓ Correctly handled missing input files"
fi

####################################################################################################

log "TEST 5: Error handling - Missing output parameter"

if "$meta_executable" \
  --files "$meta_temp_dir/sample1.bedgraph" "$meta_temp_dir/sample2.bedgraph" 2>/dev/null; then
  log "✗ Should have failed without output parameter"
  exit 1
else
  log "✓ Correctly handled missing output parameter"
fi

####################################################################################################

log "All tests completed successfully!"
