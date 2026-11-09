#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# Create test data directory
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

# --- Test Case 1: Basic sequence dictionary creation ---
log "Starting TEST 1: Basic sequence dictionary creation"

log "Generating test reference genome..."
create_test_fasta "$test_data_dir/ref.fasta" 2 500
check_file_exists "$test_data_dir/ref.fasta" "test reference genome"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/ref.fasta" \
  --output "$meta_temp_dir/ref1.dict"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/ref1.dict" "sequence dictionary file"
check_file_not_empty "$meta_temp_dir/ref1.dict" "sequence dictionary file"
check_file_matches_regex "$meta_temp_dir/ref1.dict" "^@HD" "sequence dictionary header"
check_file_matches_regex "$meta_temp_dir/ref1.dict" "^@SQ" "sequence dictionary entries"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Sequence dictionary with genome assembly and species ---
log "Starting TEST 2: Sequence dictionary with genome assembly and species"

log "Executing $meta_name with genome assembly and species options..."
"$meta_executable" \
  --input "$test_data_dir/ref.fasta" \
  --output "$meta_temp_dir/ref2.dict" \
  --genome_assembly GRCh38 \
  --species Homo_sapiens

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/ref2.dict" "sequence dictionary file with options"
check_file_not_empty "$meta_temp_dir/ref2.dict" "sequence dictionary file with options"
check_file_matches_regex "$meta_temp_dir/ref2.dict" "^@HD" "sequence dictionary header"
check_file_matches_regex "$meta_temp_dir/ref2.dict" "^@SQ" "sequence dictionary entries"
check_file_contains "$meta_temp_dir/ref2.dict" "AS:GRCh38" "genome assembly field"
check_file_contains "$meta_temp_dir/ref2.dict" "SP:Homo_sapiens" "species field"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Sequence dictionary with alternative contig names ---
log "Starting TEST 3: Sequence dictionary with alternative contig names"

log "Generating alternative names file..."
printf "seq1\tchr1\n" > "$test_data_dir/alt_names.txt"

log "Executing $meta_name with --alt_names option..."
"$meta_executable" \
  --input "$test_data_dir/ref.fasta" \
  --output "$meta_temp_dir/ref3.dict" \
  --alt_names "$test_data_dir/alt_names.txt"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/ref3.dict" "sequence dictionary file with alternative names"
check_file_not_empty "$meta_temp_dir/ref3.dict" "sequence dictionary file with alternative names"
check_file_contains "$meta_temp_dir/ref3.dict" "AN:chr1" "alternative contig name field"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Sequence dictionary limited to a number of sequences ---
log "Starting TEST 4: Sequence dictionary limited with --num_sequences"

log "Generating test reference genome with multiple sequences..."
create_test_fasta "$test_data_dir/ref_multi.fasta" 4 500
check_file_exists "$test_data_dir/ref_multi.fasta" "multi-sequence test reference genome"

log "Executing $meta_name with --num_sequences option..."
"$meta_executable" \
  --input "$test_data_dir/ref_multi.fasta" \
  --output "$meta_temp_dir/ref4.dict" \
  --num_sequences 2

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/ref4.dict" "sequence dictionary file limited by num_sequences"
# One @HD header line plus 2 @SQ entries (num_sequences)
check_file_line_count "$meta_temp_dir/ref4.dict" 3 "sequence dictionary file limited by num_sequences"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests completed successfully"
