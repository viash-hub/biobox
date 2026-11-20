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
check_file_not_exists "$test1_dir/marked_duplicates.bai" "output BAM index (should not exist without --create_index)"

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
# Check deplicate-removed BAM has fewer reads than the input BAM
input_read_count=$(count_bam_reads "$input_bam")
dedup_read_count=$(count_bam_reads "$test2_dir/dedup.bam")
if [ "$dedup_read_count" -ge "$input_read_count" ]; then
  echo "Error: Duplicate-removed BAM has $dedup_read_count reads, which is not fewer than the input BAM's $input_read_count reads." >&2
  exit 1
fi

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

# --- Test Case 5: Create an index at a different path ---
log "Starting TEST 5: --create_index with --output_index in a different directory"

test5_dir="$meta_temp_dir/test5"
test5_index_dir="$meta_temp_dir/test5_index"
mkdir -p "$test5_dir" "$test5_index_dir"

log "Executing $meta_name with --create_index and --output_index..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test5_dir/marked_duplicates.bam" \
  --metrics "$test5_dir/marked_dup_metrics.txt" \
  --create_index \
  --output_index "$test5_index_dir/marked.bai"

log "Validating TEST 5 outputs..."
check_file_exists "$test5_dir/marked_duplicates.bam" "dedup-tagged output BAM"
check_file_not_empty "$test5_dir/marked_duplicates.bam" "dedup-tagged output BAM"
check_file_exists "$test5_index_dir/marked.bai" "output BAM index"
check_file_not_empty "$test5_index_dir/marked.bai" "output BAM index"
check_file_not_exists "$test5_dir/marked_duplicates.bai" "index next to the output BAM (moved to --output_index)"
check_file_not_exists "$test5_dir/marked_duplicates.bam.bai" "index next to the output BAM (moved to --output_index)"

log "Executing $meta_name with --output_index at the path Picard writes..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test5_dir/same_path.bam" \
  --metrics "$test5_dir/same_path_metrics.txt" \
  --create_index \
  --output_index "$test5_dir/same_path.bai"

check_file_exists "$test5_dir/same_path.bai" "output BAM index (same path Picard writes)"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Index argument checks ---
log "Starting TEST 6: --create_index requires --output_index"

test6_dir="$meta_temp_dir/test6"
mkdir -p "$test6_dir"

if "$meta_executable" \
  --input "$input_bam" \
  --output "$test6_dir/fail.bam" \
  --metrics "$test6_dir/fail_metrics.txt" \
  --create_index > "$test6_dir/fail.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --create_index and no --output_index"
  exit 1
fi
check_file_contains "$test6_dir/fail.log" "Error: --create_index requires --output_index." "error message"
check_file_not_exists "$test6_dir/fail.bam" "output BAM (should not be written)"

log "Executing $meta_name with --output_index but without --create_index..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test6_dir/no_index.bam" \
  --metrics "$test6_dir/no_index_metrics.txt" \
  --output_index "$test6_dir/no_index_out.bai"

check_file_exists "$test6_dir/no_index.bam" "dedup-tagged output BAM"
check_file_not_exists "$test6_dir/no_index_out.bai" "output BAM index (ignored without --create_index)"
check_file_not_exists "$test6_dir/no_index.bai" "index next to the output BAM"

log "✅ TEST 6 completed successfully"

# --- CRAM fixtures ---
# The FASTA and its .fai sit in different directories so that the tests
# show that the component stages them together.
log "Creating reference and CRAM fixtures..."
ref_dir="$test_data_dir/reference"
ref_fai_dir="$test_data_dir/reference_fai"
mkdir -p "$ref_dir" "$ref_fai_dir"
reference="$ref_dir/reference.fasta"
create_test_reference "$reference" 1 2000
reference_fai="$ref_fai_dir/reference.fasta.fai"
mv "${reference}.fai" "$reference_fai"

input_cram="$test_data_dir/cram/input.cram"
create_test_cram "$input_bam" "$reference" "$input_cram" "$test_data_dir/cram_index/input.cram.crai"

# --- Test Case 7: CRAM input ---
log "Starting TEST 7: CRAM and BAM input with --reference_sequence"

test7_dir="$meta_temp_dir/test7"
mkdir -p "$test7_dir"

log "Executing $meta_name with a CRAM and a BAM --input..."
"$meta_executable" \
  --input "$input_cram" \
  --input "$input_bam2" \
  --output "$test7_dir/marked_duplicates.bam" \
  --metrics "$test7_dir/marked_dup_metrics.txt" \
  --reference_sequence "$reference" \
  --reference_sequence_fai "$reference_fai"

log "Validating TEST 7 outputs..."
check_file_exists "$test7_dir/marked_duplicates.bam" "dedup-tagged output BAM"
check_file_not_empty "$test7_dir/marked_duplicates.bam" "dedup-tagged output BAM"
check_file_contains "$test7_dir/marked_dup_metrics.txt" "ESTIMATED_LIBRARY_SIZE" "duplication metrics file"
expected_reads=$(( $(count_bam_reads "$input_bam") + $(count_bam_reads "$input_bam2") ))
test7_reads=$(count_bam_reads "$test7_dir/marked_duplicates.bam")
if [[ "$test7_reads" != "$expected_reads" ]]; then
  log_error "✗ Output BAM has $test7_reads reads, expected $expected_reads"
  exit 1
fi
log "✓ Output BAM has all $expected_reads input reads"

log "✅ TEST 7 completed successfully"

# --- Test Case 8: CRAM output with an index ---
log "Starting TEST 8: BAM to CRAM output with --create_index"

test8_dir="$meta_temp_dir/test8"
test8_index_dir="$meta_temp_dir/test8_index"
mkdir -p "$test8_dir" "$test8_index_dir"

log "Executing $meta_name with a CRAM --output..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test8_dir/marked_duplicates.cram" \
  --metrics "$test8_dir/marked_dup_metrics.txt" \
  --reference_sequence "$reference" \
  --reference_sequence_fai "$reference_fai" \
  --create_index \
  --output_index "$test8_index_dir/marked.cram.bai"

log "Validating TEST 8 outputs..."
check_file_exists "$test8_dir/marked_duplicates.cram" "dedup-tagged output CRAM"
check_file_not_empty "$test8_dir/marked_duplicates.cram" "dedup-tagged output CRAM"
if [[ "$(head -c 4 "$test8_dir/marked_duplicates.cram")" != "CRAM" ]]; then
  log_error "✗ Output does not have the CRAM magic bytes"
  exit 1
fi
log "✓ Output has the CRAM magic bytes"
check_file_exists "$test8_index_dir/marked.cram.bai" "output CRAM index"
check_file_not_empty "$test8_index_dir/marked.cram.bai" "output CRAM index"
if [[ "$(head -c 3 "$test8_index_dir/marked.cram.bai")" != "BAI" ]]; then
  log_error "✗ Output index does not have the BAI magic bytes"
  exit 1
fi
log "✓ Output index has the BAI magic bytes"
check_file_not_exists "$test8_dir/marked_duplicates.cram.bai" "index next to the output CRAM (moved to --output_index)"
check_file_not_exists "$test8_dir/marked_duplicates.bai" "index next to the output CRAM"
check_file_not_exists "$test8_dir/marked_duplicates.cram.crai" "CRAI index next to the output CRAM"

log "✅ TEST 8 completed successfully"

# --- Test Case 9: Reference argument checks ---
log "Starting TEST 9: CRAM requires --reference_sequence, which requires --reference_sequence_fai"

test9_dir="$meta_temp_dir/test9"
mkdir -p "$test9_dir"

log "Executing $meta_name with a CRAM --input and no --reference_sequence..."
if "$meta_executable" \
  --input "$input_cram" \
  --output "$test9_dir/cram_in.bam" \
  --metrics "$test9_dir/cram_in_metrics.txt" > "$test9_dir/cram_in.log" 2>&1; then
  log_error "✗ $meta_name did not fail with a CRAM --input and no --reference_sequence"
  exit 1
fi
check_file_contains "$test9_dir/cram_in.log" "Error: --reference_sequence is required for CRAM input or output." "error message"
check_file_not_exists "$test9_dir/cram_in.bam" "output BAM (should not be written)"

log "Executing $meta_name with a CRAM --output and no --reference_sequence..."
if "$meta_executable" \
  --input "$input_bam" \
  --output "$test9_dir/cram_out.cram" \
  --metrics "$test9_dir/cram_out_metrics.txt" > "$test9_dir/cram_out.log" 2>&1; then
  log_error "✗ $meta_name did not fail with a CRAM --output and no --reference_sequence"
  exit 1
fi
check_file_contains "$test9_dir/cram_out.log" "Error: --reference_sequence is required for CRAM input or output." "error message"
check_file_not_exists "$test9_dir/cram_out.cram" "output CRAM (should not be written)"

log "Executing $meta_name with --reference_sequence and no --reference_sequence_fai..."
if "$meta_executable" \
  --input "$input_cram" \
  --output "$test9_dir/no_fai.bam" \
  --metrics "$test9_dir/no_fai_metrics.txt" \
  --reference_sequence "$reference" > "$test9_dir/no_fai.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --reference_sequence and no --reference_sequence_fai"
  exit 1
fi
check_file_contains "$test9_dir/no_fai.log" "Error: --reference_sequence_fai must be provided when --reference_sequence is set." "error message"
check_file_not_exists "$test9_dir/no_fai.bam" "output BAM (should not be written)"

log "✅ TEST 9 completed successfully"

print_test_summary "All tests"
