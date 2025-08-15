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

# Generate test reference genome
log "Generating test reference genome..."
create_test_fasta "$test_data_dir/reference.fasta" 1 200
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"

# Build BWA index
log "Building BWA index for alignment tests..."
mkdir -p "$test_data_dir/index"
cp "$test_data_dir/reference.fasta" "$test_data_dir/index/"
bwa index "$test_data_dir/index/reference.fasta" >/dev/null 2>&1

# Verify index was created
index_files=(
  "$test_data_dir/index/reference.fasta.amb"
  "$test_data_dir/index/reference.fasta.ann"
  "$test_data_dir/index/reference.fasta.bwt"
  "$test_data_dir/index/reference.fasta.pac"
  "$test_data_dir/index/reference.fasta.sa"
)

for file in "${index_files[@]}"; do
  check_file_exists "$file" "BWA index file $(basename "$file")"
done

# Generate test FASTQ files (shorter reads for BWA aln)
log "Generating test FASTQ files for BWA aln..."
create_test_fastq "$test_data_dir/reads.fastq" 10 35
check_file_exists "$test_data_dir/reads.fastq" "test reads"

# --- Test Case 1: Basic alignment ---
log "Starting TEST 1: Basic BWA aln alignment"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --index "$test_data_dir/index/reference.fasta" \
  --reads "$test_data_dir/reads.fastq" \
  --output "$meta_temp_dir/output.sai"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output.sai" "SAI output"
check_file_not_empty "$meta_temp_dir/output.sai" "SAI output"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Custom parameters ---
log "Starting TEST 2: BWA aln with custom parameters"

log "Executing $meta_name with custom parameters..."
"$meta_executable" \
  --index "$test_data_dir/index/reference.fasta" \
  --reads "$test_data_dir/reads.fastq" \
  --output "$meta_temp_dir/custom.sai" \
  --max_diff "0.05" \
  --max_gap_opens 2

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/custom.sai" "custom SAI output"
check_file_not_empty "$meta_temp_dir/custom.sai" "custom SAI output"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Standard output ---
log "Starting TEST 3: BWA aln with stdout output"

log "Executing $meta_name with stdout output..."
stdout_output=$("$meta_executable" \
  --index "$test_data_dir/index/reference.fasta" \
  --reads "$test_data_dir/reads.fastq" 2>/dev/null)

log "Validating TEST 3 outputs..."
if [[ -n "$stdout_output" ]]; then
  log "✓ Standard output contains data"
else
  log_error "Standard output is empty"
  exit 1
fi

log "✅ TEST 3 completed successfully"

print_test_summary "All tests completed successfully"
