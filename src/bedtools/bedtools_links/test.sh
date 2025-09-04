#!/bin/bash

set -eo pipefail

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Setup test environment with strict error handling
setup_test_env

log "Starting tests for bedtools_links"

# Create test BED data
log "Creating test BED data..."
cat > "$meta_temp_dir/genes.bed" << 'EOF'
chr21	9928613	10012791	uc002yip.1	0	-
chr21	9928613	10012791	uc002yiq.1	0	-
chr21	9928613	10012791	uc002yir.1	0	-
chr21	9928613	10012791	uc010gkv.1	0	-
chr21	9928613	10061300	uc002yis.1	0	-
chr21	10042683	10120796	uc002yit.1	0	-
chr21	10042683	10120808	uc002yiu.1	0	-
chr21	10079666	10120808	uc002yiv.1	0	-
chr21	10080031	10081687	uc002yiw.1	0	-
chr21	10081660	10120796	uc002yix.2	0	-
EOF

#############################
# Test 1: Basic HTML generation
#############################
log "Starting TEST 1: Basic HTML generation"

log "Executing bedtools_links with default parameters..."
"$meta_executable" \
  --input "$meta_temp_dir/genes.bed" \
  --output "$meta_temp_dir/output1.html"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output1.html" "HTML output"
check_file_not_empty "$meta_temp_dir/output1.html" "HTML output"
check_file_contains "$meta_temp_dir/output1.html" "uc002yip.1" "HTML output"

log "✅ TEST 1 completed successfully"

#############################
# Test 2: Custom base URL
#############################
log "Starting TEST 2: Custom base URL"

log "Executing bedtools_links with custom base URL..."
"$meta_executable" \
  --input "$meta_temp_dir/genes.bed" \
  --output "$meta_temp_dir/output2.html" \
  --base_url "http://genome.ucsc.edu"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output2.html" "HTML output with custom base URL"
check_file_not_empty "$meta_temp_dir/output2.html" "HTML output with custom base URL"
check_file_contains "$meta_temp_dir/output2.html" "uc002yip.1" "HTML output with custom base URL"

log "✅ TEST 2 completed successfully"

#############################
# Test 3: Custom organism and database
#############################
log "Starting TEST 3: Custom organism and database"

log "Executing bedtools_links with custom organism and database..."
"$meta_executable" \
  --input "$meta_temp_dir/genes.bed" \
  --output "$meta_temp_dir/output3.html" \
  --base_url "http://genome.ucsc.edu" \
  --organism "mouse" \
  --database "mm9"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output3.html" "HTML output with mouse/mm9"
check_file_not_empty "$meta_temp_dir/output3.html" "HTML output with mouse/mm9"
check_file_contains "$meta_temp_dir/output3.html" "uc002yip.1" "HTML output with mouse/mm9"

log "✅ TEST 3 completed successfully"

log "All tests completed successfully!"
