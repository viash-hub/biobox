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
create_test_fasta "$test_data_dir/reference.fasta" 1 500
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

# Generate test FASTQ files
log "Generating test FASTQ files..."
create_test_fastq "$test_data_dir/reads_single.fastq" 15 60
create_test_fastq "$test_data_dir/reads_R1.fastq" 15 60
create_test_fastq "$test_data_dir/reads_R2.fastq" 15 60
check_file_exists "$test_data_dir/reads_single.fastq" "single-end reads"
check_file_exists "$test_data_dir/reads_R1.fastq" "paired-end R1 reads"
check_file_exists "$test_data_dir/reads_R2.fastq" "paired-end R2 reads"

# --- Test Case 1: Single-end alignment ---
log "Starting TEST 1: Single-end BWA MEM alignment"

log "Executing $meta_name with single-end reads..."
"$meta_executable" \
  --index "$test_data_dir/index/reference.fasta" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/single_end.sam"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/single_end.sam" "single-end SAM output"
check_file_not_empty "$meta_temp_dir/single_end.sam" "single-end SAM output"

# Check SAM format headers
if head -5 "$meta_temp_dir/single_end.sam" | grep -q "^@"; then
  log "✓ SAM file contains proper headers"
else
  log_error "SAM file does not contain proper headers"
  exit 1
fi

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Paired-end alignment ---
log "Starting TEST 2: Paired-end BWA MEM alignment"

log "Executing $meta_name with paired-end reads..."
"$meta_executable" \
  --index "$test_data_dir/index/reference.fasta" \
  --reads1 "$test_data_dir/reads_R1.fastq" \
  --reads2 "$test_data_dir/reads_R2.fastq" \
  --output "$meta_temp_dir/paired_end.sam"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/paired_end.sam" "paired-end SAM output"
check_file_not_empty "$meta_temp_dir/paired_end.sam" "paired-end SAM output"

# Check SAM format headers
if head -5 "$meta_temp_dir/paired_end.sam" | grep -q "^@"; then
  log "✓ SAM file contains proper headers"
else
  log_error "SAM file does not contain proper headers"
  exit 1
fi

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Advanced parameters ---
log "Starting TEST 3: BWA MEM with advanced parameters"

log "Executing $meta_name with advanced parameters..."
"$meta_executable" \
  --index "$test_data_dir/index/reference.fasta" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/advanced.sam" \
  --threads 2 \
  --min_seed_length 15

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/advanced.sam" "advanced SAM output"
check_file_not_empty "$meta_temp_dir/advanced.sam" "advanced SAM output"

log "✅ TEST 3 completed successfully"

print_test_summary "All tests completed successfully"
