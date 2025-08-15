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

# Generate test genome for indexing
log "Generating test genome..."
create_test_fasta "$test_data_dir/genome.fasta" 2 1000
check_file_exists "$test_data_dir/genome.fasta" "test genome"

# Generate test FASTQ files
log "Generating test FASTQ files..."
create_test_fastq "$test_data_dir/reads_single.fastq" 20 50
create_test_fastq "$test_data_dir/reads_R1.fastq" 20 50
create_test_fastq "$test_data_dir/reads_R2.fastq" 20 50
check_file_exists "$test_data_dir/reads_single.fastq" "single-end reads"
check_file_exists "$test_data_dir/reads_R1.fastq" "paired-end R1 reads"
check_file_exists "$test_data_dir/reads_R2.fastq" "paired-end R2 reads"

# Build bowtie2 index
log "Building bowtie2 index..."
mkdir -p "$test_data_dir/index"
bowtie2-build "$test_data_dir/genome.fasta" "$test_data_dir/index/genome" >/dev/null 2>&1
check_file_exists "$test_data_dir/index/genome.1.bt2" "bowtie2 index file"

# --- Test Case 1: Single-end alignment ---
log "Starting TEST 1: Single-end alignment"

log "Executing $meta_name with single-end reads..."
"$meta_executable" \
  --index "$test_data_dir/index/genome" \
  --unpaired "$test_data_dir/reads_single.fastq" \
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
log "Starting TEST 2: Paired-end alignment"

log "Executing $meta_name with paired-end reads..."
"$meta_executable" \
  --index "$test_data_dir/index/genome" \
  --mate1 "$test_data_dir/reads_R1.fastq" \
  --mate2 "$test_data_dir/reads_R2.fastq" \
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

# --- Test Case 3: Advanced alignment parameters ---
log "Starting TEST 3: Advanced alignment parameters"

log "Executing $meta_name with advanced parameters..."
"$meta_executable" \
  --index "$test_data_dir/index/genome" \
  --unpaired "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/advanced.sam" \
  --threads 2 \
  --very_sensitive \
  --no_unal

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/advanced.sam" "advanced SAM output"
check_file_not_empty "$meta_temp_dir/advanced.sam" "advanced SAM output"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Output format options ---
log "Starting TEST 4: Output format options"

log "Executing $meta_name with BAM output..."
"$meta_executable" \
  --index "$test_data_dir/index/genome" \
  --unpaired "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/output.bam" \
  --threads 2

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output.bam" "BAM output"
check_file_not_empty "$meta_temp_dir/output.bam" "BAM output"

# Verify BAM format by checking if samtools can read it
if command -v samtools >/dev/null 2>&1; then
  if samtools view -H "$meta_temp_dir/output.bam" >/dev/null 2>&1; then
    log "✓ BAM file format is valid"
  else
    log_error "BAM file format is invalid"
    exit 1
  fi
else
  log_warn "samtools not available, skipping BAM format validation"
fi

log "✅ TEST 4 completed successfully"

print_test_summary "All tests completed successfully"
