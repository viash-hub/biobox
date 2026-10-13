#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers and the GATK4-specific helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

# --- Shared fixtures ---
# First BAM: a duplicate pair (read1/read2 at position 1) and two reads at
# distinct positions (500, 800)
create_test_sam "$test_data_dir/input.sam" seq1 2000 rg1 sample1 1 1 500 800
input_bam="$test_data_dir/input.bam"
sort_and_index_bam "$test_data_dir/input.sam" "$input_bam" coordinate
check_file_not_empty "$input_bam" "coordinate-sorted input BAM"

# Second BAM: a duplicate of the first BAM's position 1 pair and one read at
# a new position (1300)
create_test_sam "$test_data_dir/input2.sam" seq1 2000 rg1 sample1 1 1300
input_bam2="$test_data_dir/input2.bam"
sort_and_index_bam "$test_data_dir/input2.sam" "$input_bam2" coordinate
check_file_not_empty "$input_bam2" "second coordinate-sorted input BAM"

# Reference and CRAM. The FASTA and its .fai sit in different directories so
# that the tests show that the component stages them together.
ref_dir="$test_data_dir/reference"
ref_fai_dir="$test_data_dir/reference_fai"
mkdir -p "$ref_dir" "$ref_fai_dir"
reference="$ref_dir/reference.fasta"
create_test_reference "$reference" 1 2000
reference_fai="$ref_fai_dir/reference.fasta.fai"
mv "${reference}.fai" "$reference_fai"

input_cram="$test_data_dir/cram/input.cram"
create_test_cram "$input_bam" "$reference" "$input_cram" "$test_data_dir/cram_index/input.cram.crai"

input_reads=$(count_bam_reads "$input_bam")
input2_reads=$(count_bam_reads "$input_bam2")

# --- Test Case 1: Basic run (defaults) ---
log "Starting TEST 1: Basic run (defaults)"

test1_dir="$meta_temp_dir/test1"
mkdir -p "$test1_dir"

# --output_index is set without --create_index, it must be ignored
"$meta_executable" \
  --input "$input_bam" \
  --output "$test1_dir/marked_duplicates.bam" \
  --metrics "$test1_dir/marked_dup_metrics.txt" \
  --output_index "$test1_dir/marked_out.bai"

log "Validating TEST 1 outputs..."
check_file_not_empty "$test1_dir/marked_duplicates.bam" "dedup-tagged output BAM"
check_file_not_empty "$test1_dir/marked_dup_metrics.txt" "duplication metrics file"
check_file_contains "$test1_dir/marked_dup_metrics.txt" "ESTIMATED_LIBRARY_SIZE" "duplication metrics file"
check_file_not_exists "$test1_dir/marked_duplicates.bai" "output BAM index (should not exist without --create_index)"
check_file_not_exists "$test1_dir/marked_out.bai" "output BAM index (--output_index ignored without --create_index)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs ---
log "Starting TEST 2: CRAM and BAM input with --reference_sequence"

test2_dir="$meta_temp_dir/test2"
mkdir -p "$test2_dir"

"$meta_executable" \
  --input "$input_cram" \
  --input "$input_bam2" \
  --output "$test2_dir/marked_duplicates.bam" \
  --metrics "$test2_dir/marked_dup_metrics.txt" \
  --reference_sequence "$reference" \
  --reference_sequence_fai "$reference_fai"

log "Validating TEST 2 outputs..."
check_file_not_empty "$test2_dir/marked_duplicates.bam" "dedup-tagged output BAM"
check_file_contains "$test2_dir/marked_dup_metrics.txt" "ESTIMATED_LIBRARY_SIZE" "duplication metrics file"
expected_reads=$((input_reads + input2_reads))
test2_reads=$(count_bam_reads "$test2_dir/marked_duplicates.bam")
if [[ "$test2_reads" != "$expected_reads" ]]; then
  log_error "✗ Output BAM has $test2_reads reads, expected $expected_reads"
  exit 1
fi
log "✓ Output BAM has all $expected_reads input reads"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Outputs and index options ---
log "Starting TEST 3: Outputs and index options"

test3_dir="$meta_temp_dir/test3"
test3_index_dir="$meta_temp_dir/test3_index"
mkdir -p "$test3_dir" "$test3_index_dir"

# Run A: index in another directory, with the tagging and metadata options.
# The options have a visible effect in the output BAM, except
# --tag_duplicate_set_members.
log "Executing $meta_name with --create_index, tagging and metadata options..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test3_dir/tagged.bam" \
  --metrics "$test3_dir/tagged_metrics.txt" \
  --create_index \
  --output_index "$test3_index_dir/tagged.bai" \
  --tagging_policy All \
  --tag_duplicate_set_members \
  --comment "first comment" \
  --comment "second comment" \
  --program_record_id MarkDuplicatesTest \
  --program_group_name MarkDuplicatesTestGroup

check_file_not_empty "$test3_dir/tagged.bam" "tagged output BAM"
check_file_contains "$test3_dir/tagged_metrics.txt" "ESTIMATED_LIBRARY_SIZE" "duplication metrics file"
check_index_moved "$test3_index_dir/tagged.bai" "$test3_dir/tagged.bai" "output BAM index"
check_file_not_empty "$test3_index_dir/tagged.bai" "output BAM index"
check_file_not_exists "$test3_dir/tagged.bam.bai" "index next to the output BAM (moved to --output_index)"

samtools view -H "$test3_dir/tagged.bam" > "$test3_dir/tagged_header.sam"
samtools view "$test3_dir/tagged.bam" > "$test3_dir/tagged_reads.sam"
check_file_contains "$test3_dir/tagged_header.sam" $'@CO\tfirst comment' "header comment (--comment)"
check_file_contains "$test3_dir/tagged_header.sam" $'@CO\tsecond comment' "header comment (--comment, second value)"
check_file_contains "$test3_dir/tagged_header.sam" $'@PG\tID:MarkDuplicatesTest\t' "program record ID (--program_record_id)"
check_file_contains "$test3_dir/tagged_header.sam" "PN:MarkDuplicatesTestGroup" "program group name (--program_group_name)"
check_file_contains "$test3_dir/tagged_reads.sam" "DT:Z:LB" "duplicate type tag (--tagging_policy All)"
# Picard only adds the duplicate set tags (DS, DI) to read pairs, and these
# reads are single-end, so check the Picard command line in the header
check_file_contains "$test3_dir/tagged_header.sam" "TAG_DUPLICATE_SET_MEMBERS true" "@PG command line (--tag_duplicate_set_members)"

# Run B: --remove_duplicates, index at the path Picard writes
log "Executing $meta_name with --remove_duplicates and --output_index at the path Picard writes..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test3_dir/dedup.bam" \
  --metrics "$test3_dir/dedup_metrics.txt" \
  --remove_duplicates \
  --create_index \
  --output_index "$test3_dir/dedup.bai"

check_file_not_empty "$test3_dir/dedup.bam" "duplicate-removed output BAM"
check_file_contains "$test3_dir/dedup_metrics.txt" "ESTIMATED_LIBRARY_SIZE" "duplication metrics file"
check_file_not_empty "$test3_dir/dedup.bai" "output BAM index (same path Picard writes)"
dedup_reads=$(count_bam_reads "$test3_dir/dedup.bam")
if [[ "$dedup_reads" -ge "$input_reads" ]]; then
  log_error "✗ Duplicate-removed BAM has $dedup_reads reads, which is not fewer than the input BAM's $input_reads reads"
  exit 1
fi
log "✓ Duplicate-removed BAM has fewer reads ($dedup_reads) than the input BAM ($input_reads)"

# Run C: CRAM output with an index in another directory
log "Executing $meta_name with a CRAM --output and --create_index..."
"$meta_executable" \
  --input "$input_bam" \
  --output "$test3_dir/marked_duplicates.cram" \
  --metrics "$test3_dir/cram_metrics.txt" \
  --reference_sequence "$reference" \
  --reference_sequence_fai "$reference_fai" \
  --create_index \
  --output_index "$test3_index_dir/marked.cram.bai"

check_file_not_empty "$test3_dir/marked_duplicates.cram" "dedup-tagged output CRAM"
if [[ "$(head -c 4 "$test3_dir/marked_duplicates.cram")" != "CRAM" ]]; then
  log_error "✗ Output does not have the CRAM magic bytes"
  exit 1
fi
log "✓ Output has the CRAM magic bytes"
check_index_moved "$test3_index_dir/marked.cram.bai" "$test3_dir/marked_duplicates.cram.bai" "output CRAM index"
check_file_not_empty "$test3_index_dir/marked.cram.bai" "output CRAM index"
if [[ "$(head -c 3 "$test3_index_dir/marked.cram.bai")" != "BAI" ]]; then
  log_error "✗ Output index does not have the BAI magic bytes"
  exit 1
fi
log "✓ Output index has the BAI magic bytes"
check_file_not_exists "$test3_dir/marked_duplicates.bai" "index next to the output CRAM"
check_file_not_exists "$test3_dir/marked_duplicates.cram.crai" "CRAI index next to the output CRAM"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Argument errors ---
log "Starting TEST 4: Argument errors"

test4_dir="$meta_temp_dir/test4"
mkdir -p "$test4_dir"

log "Executing $meta_name with --create_index and no --output_index..."
check_fails_with "Error: --create_index requires --output_index." "$test4_dir/no_output_index.log" \
  "$meta_executable" \
  --input "$input_bam" \
  --output "$test4_dir/no_output_index.bam" \
  --metrics "$test4_dir/no_output_index_metrics.txt" \
  --create_index
check_file_not_exists "$test4_dir/no_output_index.bam" "output BAM (should not be written)"

log "Executing $meta_name with a CRAM --input and no --reference_sequence..."
check_fails_with "Error: --reference_sequence is required for CRAM input or output." "$test4_dir/cram_in.log" \
  "$meta_executable" \
  --input "$input_cram" \
  --output "$test4_dir/cram_in.bam" \
  --metrics "$test4_dir/cram_in_metrics.txt"
check_file_not_exists "$test4_dir/cram_in.bam" "output BAM (should not be written)"

log "Executing $meta_name with a CRAM --output and no --reference_sequence..."
check_fails_with "Error: --reference_sequence is required for CRAM input or output." "$test4_dir/cram_out.log" \
  "$meta_executable" \
  --input "$input_bam" \
  --output "$test4_dir/cram_out.cram" \
  --metrics "$test4_dir/cram_out_metrics.txt"
check_file_not_exists "$test4_dir/cram_out.cram" "output CRAM (should not be written)"

log "Executing $meta_name with --reference_sequence and no --reference_sequence_fai..."
check_fails_with "Error: --reference_sequence_fai must be provided when --reference_sequence is set." "$test4_dir/no_fai.log" \
  "$meta_executable" \
  --input "$input_cram" \
  --output "$test4_dir/no_fai.bam" \
  --metrics "$test4_dir/no_fai_metrics.txt" \
  --reference_sequence "$reference"
check_file_not_exists "$test4_dir/no_fai.bam" "output BAM (should not be written)"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests"
