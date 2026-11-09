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

# Build a standard minibwa index
log "Building minibwa index for alignment tests..."
mkdir -p "$test_data_dir/index"
minibwa index "$test_data_dir/reference.fasta" "$test_data_dir/index/reference" >/dev/null 2>&1

index_files=(
  "$test_data_dir/index/reference.l2b"
  "$test_data_dir/index/reference.mbw"
)

for file in "${index_files[@]}"; do
  check_file_exists "$file" "minibwa index file $(basename "$file")"
done

# Build a bisulfite minibwa index
log "Building bisulfite minibwa index for bisulfite alignment test..."
mkdir -p "$test_data_dir/bisulfite_index"
minibwa index --meth "$test_data_dir/reference.fasta" "$test_data_dir/bisulfite_index/reference" >/dev/null 2>&1
check_file_exists "$test_data_dir/bisulfite_index/reference.meth.mbw" "bisulfite minibwa index file"

# Generate test FASTQ files
log "Generating test FASTQ files..."
create_test_fastq "$test_data_dir/reads_single.fastq" 15 60
create_test_fastq "$test_data_dir/reads_R1.fastq" 15 60
create_test_fastq "$test_data_dir/reads_R2.fastq" 15 60
check_file_exists "$test_data_dir/reads_single.fastq" "single-end reads"
check_file_exists "$test_data_dir/reads_R1.fastq" "paired-end R1 reads"
check_file_exists "$test_data_dir/reads_R2.fastq" "paired-end R2 reads"

# --- Test Case 1: Single-end alignment ---
log "Starting TEST 1: Single-end minibwa map alignment"

log "Executing $meta_name with single-end reads..."
"$meta_executable" \
  --index "$test_data_dir/index" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/single_end.sam"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/single_end.sam" "single-end SAM output"
check_file_not_empty "$meta_temp_dir/single_end.sam" "single-end SAM output"

if head -5 "$meta_temp_dir/single_end.sam" | grep -q "^@"; then
  log "✓ SAM file contains proper headers"
else
  log_error "SAM file does not contain proper headers"
  exit 1
fi

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Paired-end alignment ---
log "Starting TEST 2: Paired-end minibwa map alignment"

log "Executing $meta_name with paired-end reads..."
"$meta_executable" \
  --index "$test_data_dir/index" \
  --reads1 "$test_data_dir/reads_R1.fastq" \
  --reads2 "$test_data_dir/reads_R2.fastq" \
  --output "$meta_temp_dir/paired_end.sam"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/paired_end.sam" "paired-end SAM output"
check_file_not_empty "$meta_temp_dir/paired_end.sam" "paired-end SAM output"

if head -5 "$meta_temp_dir/paired_end.sam" | grep -q "^@"; then
  log "✓ SAM file contains proper headers"
else
  log_error "SAM file does not contain proper headers"
  exit 1
fi

log "✅ TEST 2 completed successfully"

# --- Test Case 3: PAF output ---
log "Starting TEST 3: minibwa map with --paf_output"

log "Executing $meta_name with --paf_output..."
"$meta_executable" \
  --index "$test_data_dir/index" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/output.paf" \
  --paf_output

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output.paf" "PAF output"
check_file_not_empty "$meta_temp_dir/output.paf" "PAF output"

if head -1 "$meta_temp_dir/output.paf" | grep -q "^@"; then
  log_error "PAF output should not contain SAM-style headers"
  exit 1
else
  log "✓ PAF output does not contain SAM-style headers"
fi

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Bisulfite alignment ---
log "Starting TEST 4: minibwa map with --meth"

log "Executing $meta_name with --meth..."
"$meta_executable" \
  --index "$test_data_dir/bisulfite_index" \
  --reads1 "$test_data_dir/reads_R1.fastq" \
  --reads2 "$test_data_dir/reads_R2.fastq" \
  --meth \
  --output "$meta_temp_dir/bisulfite.sam"

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/bisulfite.sam" "bisulfite SAM output"
check_file_not_empty "$meta_temp_dir/bisulfite.sam" "bisulfite SAM output"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests"
