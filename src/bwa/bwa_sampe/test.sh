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

# Build a BWA index in a dedicated directory, laid out exactly like the output of
# `bwa_index`: only the index files, no reference FASTA next to them. This is the
# layout that has to work, because a directory is what gets mounted (Docker) or
# staged (Nextflow) in full.
log "Building BWA index for sampe tests..."
mkdir -p "$test_data_dir/index"
bwa index -p "$test_data_dir/index/reference" "$test_data_dir/reference.fasta" >/dev/null 2>&1

index_files=(
  "$test_data_dir/index/reference.amb"
  "$test_data_dir/index/reference.ann"
  "$test_data_dir/index/reference.bwt"
  "$test_data_dir/index/reference.pac"
  "$test_data_dir/index/reference.sa"
)

for file in "${index_files[@]}"; do
  check_file_exists "$file" "BWA index file $(basename "$file")"
done

check_file_not_exists "$test_data_dir/index/reference.fasta" "reference FASTA in the index directory"

# Generate test FASTQ files (shorter reads for BWA aln)
log "Generating test FASTQ files for BWA aln..."
create_test_fastq "$test_data_dir/reads_R1.fastq" 8 35
create_test_fastq "$test_data_dir/reads_R2.fastq" 8 35
check_file_exists "$test_data_dir/reads_R1.fastq" "R1 reads"
check_file_exists "$test_data_dir/reads_R2.fastq" "R2 reads"

# Generate SAI files using BWA aln
log "Generating SAI files for sampe test..."
bwa aln "$test_data_dir/index/reference" "$test_data_dir/reads_R1.fastq" > "$test_data_dir/reads_R1.sai" 2>/dev/null
bwa aln "$test_data_dir/index/reference" "$test_data_dir/reads_R2.fastq" > "$test_data_dir/reads_R2.sai" 2>/dev/null

check_file_exists "$test_data_dir/reads_R1.sai" "R1 SAI file"
check_file_exists "$test_data_dir/reads_R2.sai" "R2 SAI file"
check_file_not_empty "$test_data_dir/reads_R1.sai" "R1 SAI file"
check_file_not_empty "$test_data_dir/reads_R2.sai" "R2 SAI file"

# --- Test Case 1: Basic paired-end SAM generation ---
log "Starting TEST 1: Basic BWA sampe"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --index "$test_data_dir/index" \
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
  --index "$test_data_dir/index" \
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

# --- Test Case 3: Index prefix unrelated to the reference file name ---
log "Starting TEST 3: BWA sampe with a custom index prefix"

# The prefix must be discovered from the directory contents, not derived from the
# name of the directory or of the reference FASTA.
log "Building BWA index with a custom prefix..."
mkdir -p "$test_data_dir/custom_index"
bwa index -p "$test_data_dir/custom_index/my_genome_v2" "$test_data_dir/reference.fasta" >/dev/null 2>&1
check_file_exists "$test_data_dir/custom_index/my_genome_v2.bwt" "custom-prefixed index file"

log "Executing $meta_name with a custom index prefix..."
"$meta_executable" \
  --index "$test_data_dir/custom_index" \
  --sai1 "$test_data_dir/reads_R1.sai" \
  --sai2 "$test_data_dir/reads_R2.sai" \
  --reads1 "$test_data_dir/reads_R1.fastq" \
  --reads2 "$test_data_dir/reads_R2.fastq" \
  --output "$meta_temp_dir/custom_prefix.sam"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/custom_prefix.sam" "custom-prefix SAM output"
check_file_not_empty "$meta_temp_dir/custom_prefix.sam" "custom-prefix SAM output"
check_file_contains "$meta_temp_dir/custom_prefix.sam" "@SQ" "custom-prefix SAM output"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Index directory without an index ---
log "Starting TEST 4: BWA sampe with an index directory that holds no index"

mkdir -p "$test_data_dir/empty_index"
cp "$test_data_dir/reference.fasta" "$test_data_dir/empty_index/"

log "Executing $meta_name with an index directory that holds no index..."
if "$meta_executable" \
  --index "$test_data_dir/empty_index" \
  --sai1 "$test_data_dir/reads_R1.sai" \
  --sai2 "$test_data_dir/reads_R2.sai" \
  --reads1 "$test_data_dir/reads_R1.fastq" \
  --reads2 "$test_data_dir/reads_R2.fastq" \
  --output "$meta_temp_dir/should_not_exist.sam" 2> "$meta_temp_dir/no_index.err"; then
  log_error "Expected a non-zero exit code for an index directory without an index"
  exit 1
fi

log "Validating TEST 4 outputs..."
check_file_contains "$meta_temp_dir/no_index.err" "No BWA index" "error message"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests completed successfully"
