#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers and the GATK4-specific helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# Create test data directory
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

# Generate a synthetic SAM fixture with a duplicate read pair: read1/read2
# share the same position (1), read3/read4 sit at distinct positions.
sam_file="$test_data_dir/input.sam"
create_test_sam "$sam_file" seq1 2000 rg1 sample1 1 1 500 800
check_file_exists "$sam_file" "synthetic SAM fixture"

log "Sorting and indexing synthetic SAM fixture..."
input_bam="$test_data_dir/input.bam"
sort_and_index_bam "$sam_file" "$input_bam" coordinate
check_file_exists "$input_bam" "coordinate-sorted input BAM"
check_file_not_empty "$input_bam" "coordinate-sorted input BAM"

# --- Test Case 1: Basic duplicate marking ---
log "Starting TEST 1: Basic duplicate marking"

test1_dir="$meta_temp_dir/test1"
mkdir -p "$test1_dir"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test1_dir/marked_duplicates.bam" \
  --metrics "$test1_dir/marked_dup_metrics.txt"

log "Validating TEST 1 outputs..."
check_file_exists "$test1_dir/marked_duplicates.bam" "dedup-tagged output BAM"
check_file_not_empty "$test1_dir/marked_duplicates.bam" "dedup-tagged output BAM"
check_file_exists "$test1_dir/marked_dup_metrics.txt" "duplication metrics file"
check_file_not_empty "$test1_dir/marked_dup_metrics.txt" "duplication metrics file"
check_file_contains "$test1_dir/marked_dup_metrics.txt" "ESTIMATED_LIBRARY_SIZE" "duplication metrics file"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Remove duplicates ---
log "Starting TEST 2: Remove duplicates"

test2_dir="$meta_temp_dir/test2"
mkdir -p "$test2_dir"

log "Executing $meta_name with --remove_duplicates..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test2_dir/dedup.bam" \
  --metrics "$test2_dir/dedup_metrics.txt" \
  --remove_duplicates

log "Validating TEST 2 outputs..."
check_file_exists "$test2_dir/dedup.bam" "duplicate-removed output BAM"
check_file_not_empty "$test2_dir/dedup.bam" "duplicate-removed output BAM"
check_file_exists "$test2_dir/dedup_metrics.txt" "duplication metrics file"
check_file_not_empty "$test2_dir/dedup_metrics.txt" "duplication metrics file"
check_file_contains "$test2_dir/dedup_metrics.txt" "ESTIMATED_LIBRARY_SIZE" "duplication metrics file"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Tagging policy, duplicate set tagging and metadata options ---
log "Starting TEST 3: Tagging policy, duplicate set tagging and metadata options"

test3_dir="$meta_temp_dir/test3"
mkdir -p "$test3_dir"

log "Executing $meta_name with --tagging_policy, --tag_duplicate_set_members, --comment and program record options..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test3_dir/marked_duplicates.bam" \
  --metrics "$test3_dir/marked_dup_metrics.txt" \
  --tagging_policy All \
  --tag_duplicate_set_members \
  --comment "first comment" \
  --comment "second comment" \
  --program_record_id MarkDuplicatesTest \
  --program_group_name MarkDuplicatesTest

log "Validating TEST 3 outputs..."
check_file_exists "$test3_dir/marked_duplicates.bam" "tagged output BAM"
check_file_not_empty "$test3_dir/marked_duplicates.bam" "tagged output BAM"
check_file_exists "$test3_dir/marked_dup_metrics.txt" "duplication metrics file"
check_file_not_empty "$test3_dir/marked_dup_metrics.txt" "duplication metrics file"
check_file_contains "$test3_dir/marked_dup_metrics.txt" "ESTIMATED_LIBRARY_SIZE" "duplication metrics file"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Multiple --input files (merge and de-duplicate together) ---
log "Starting TEST 4: Multiple --input files"

# Build a second BAM sharing a duplicate position (1) with the first file's
# read1/read2 pair, plus one read at a new position
sam_file2="$test_data_dir/input2.sam"
create_test_sam "$sam_file2" seq1 2000 rg1 sample1 1 1300
check_file_exists "$sam_file2" "second synthetic SAM fixture"

input_bam2="$test_data_dir/input2.bam"
sort_and_index_bam "$sam_file2" "$input_bam2" coordinate
check_file_exists "$input_bam2" "second coordinate-sorted input BAM"

test4_dir="$meta_temp_dir/test4"
mkdir -p "$test4_dir"

log "Executing $meta_name with two --input files..."
"$meta_executable" \
  --input "$input_bam" \
  --input "$input_bam2" \
  --output "$test4_dir/marked_duplicates.bam" \
  --metrics "$test4_dir/marked_dup_metrics.txt"

log "Validating TEST 4 outputs..."
check_file_exists "$test4_dir/marked_duplicates.bam" "merged dedup-tagged output BAM"
check_file_not_empty "$test4_dir/marked_duplicates.bam" "merged dedup-tagged output BAM"
check_file_exists "$test4_dir/marked_dup_metrics.txt" "merged duplication metrics file"
check_file_not_empty "$test4_dir/marked_dup_metrics.txt" "merged duplication metrics file"
check_file_contains "$test4_dir/marked_dup_metrics.txt" "ESTIMATED_LIBRARY_SIZE" "merged duplication metrics file"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests"
