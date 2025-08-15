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
log "Building BWA index for sampe tests..."
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
create_test_fastq "$test_data_dir/reads_R1.fastq" 8 35
create_test_fastq "$test_data_dir/reads_R2.fastq" 8 35
check_file_exists "$test_data_dir/reads_R1.fastq" "R1 reads"
check_file_exists "$test_data_dir/reads_R2.fastq" "R2 reads"

# Generate SAI files using BWA aln
log "Generating SAI files for sampe test..."
bwa aln "$test_data_dir/index/reference.fasta" "$test_data_dir/reads_R1.fastq" > "$test_data_dir/reads_R1.sai" 2>/dev/null
bwa aln "$test_data_dir/index/reference.fasta" "$test_data_dir/reads_R2.fastq" > "$test_data_dir/reads_R2.sai" 2>/dev/null

check_file_exists "$test_data_dir/reads_R1.sai" "R1 SAI file"
check_file_exists "$test_data_dir/reads_R2.sai" "R2 SAI file"
check_file_not_empty "$test_data_dir/reads_R1.sai" "R1 SAI file"
check_file_not_empty "$test_data_dir/reads_R2.sai" "R2 SAI file"

# --- Test Case 1: Basic paired-end SAM generation ---
log "Starting TEST 1: Basic BWA sampe"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --index "$test_data_dir/index/reference.fasta" \
  --sai1 "$test_data_dir/reads_R1.sai" \
  --sai2 "$test_data_dir/reads_R2.sai" \
  --reads1 "$test_data_dir/reads_R1.fastq" \
  --reads2 "$test_data_dir/reads_R2.fastq" \
  --output "$meta_temp_dir/paired_end.sam"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/paired_end.sam" "paired-end SAM output"
check_file_not_empty "$meta_temp_dir/paired_end.sam" "paired-end SAM output"

# Check SAM format headers
if head -5 "$meta_temp_dir/paired_end.sam" | grep -q "^@"; then
  log "✓ SAM file contains proper headers"
else
  log_error "SAM file does not contain proper headers"
  exit 1
fi

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Custom parameters ---
log "Starting TEST 2: BWA sampe with custom parameters"

log "Executing $meta_name with custom parameters..."
"$meta_executable" \
  --index "$test_data_dir/index/reference.fasta" \
  --sai1 "$test_data_dir/reads_R1.sai" \
  --sai2 "$test_data_dir/reads_R2.sai" \
  --reads1 "$test_data_dir/reads_R1.fastq" \
  --reads2 "$test_data_dir/reads_R2.fastq" \
  --output "$meta_temp_dir/custom.sam" \
  --max_insert_size 800

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/custom.sam" "custom SAM output"
check_file_not_empty "$meta_temp_dir/custom.sam" "custom SAM output"

# Check SAM format headers
if head -5 "$meta_temp_dir/custom.sam" | grep -q "^@"; then
  log "✓ SAM file contains proper headers"
else
  log_error "SAM file does not contain proper headers"
  exit 1
fi

log "✅ TEST 2 completed successfully"

print_test_summary "All tests completed successfully"
