#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

# --- Shared fixtures ---
log "Generating test reference genome with 4 sequences..."
create_test_fasta "$test_data_dir/ref.fasta" 4 500
check_file_exists "$test_data_dir/ref.fasta" "test reference genome"

printf "seq1\tchr1\n" > "$test_data_dir/alt_names.txt"

# --- Test Case 1: Basic run (defaults) ---
log "Starting TEST 1: Basic run (defaults)"

"$meta_executable" \
  --input "$test_data_dir/ref.fasta" \
  --output "$meta_temp_dir/basic.dict"

log "Validating TEST 1 outputs..."
check_file_not_empty "$meta_temp_dir/basic.dict" "sequence dictionary file"
check_file_matches_regex "$meta_temp_dir/basic.dict" "^@HD" "sequence dictionary header"
check_file_matches_regex "$meta_temp_dir/basic.dict" "^@SQ" "sequence dictionary entries"
# One @HD header line plus 4 @SQ entries
check_file_line_count "$meta_temp_dir/basic.dict" 5 "sequence dictionary file"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Tool options ---
log "Starting TEST 2: Tool options"
"$meta_executable" \
  --input "$test_data_dir/ref.fasta" \
  --output "$meta_temp_dir/options.dict" \
  --genome_assembly GRCh38 \
  --species Homo_sapiens \
  --alt_names "$test_data_dir/alt_names.txt" \
  --num_sequences 2

log "Validating TEST 2 outputs..."
check_file_matches_regex "$meta_temp_dir/options.dict" "^@HD" "sequence dictionary header"
check_file_contains "$meta_temp_dir/options.dict" "AS:GRCh38" "genome assembly field"
check_file_contains "$meta_temp_dir/options.dict" "SP:Homo_sapiens" "species field"
check_file_contains "$meta_temp_dir/options.dict" "AN:chr1" "alternative contig name field"
# One @HD header line plus 2 @SQ entries (--num_sequences 2)
check_file_line_count "$meta_temp_dir/options.dict" 3 "sequence dictionary file limited by --num_sequences"

log "✅ TEST 2 completed successfully"

print_test_summary "All tests"
