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

# Build a BWA-MEM2 index in a dedicated directory, laid out exactly like the
# output of `bwa_mem2_index`: only the index files, no reference FASTA next to
# them. This is the layout that has to work, because a directory is what gets
# mounted (Docker) or staged (Nextflow) in full.
log "Building BWA index for alignment tests..."
mkdir -p "$test_data_dir/index"
bwa-mem2 index -p "$test_data_dir/index/reference" "$test_data_dir/reference.fasta" >/dev/null 2>&1

index_files=(
  "$test_data_dir/index/reference.0123"
  "$test_data_dir/index/reference.amb"
  "$test_data_dir/index/reference.ann"
  "$test_data_dir/index/reference.bwt.2bit.64"
  "$test_data_dir/index/reference.pac"
)

for file in "${index_files[@]}"; do
  check_file_exists "$file" "BWA index file $(basename "$file")"
done

check_file_not_exists "$test_data_dir/index/reference.fasta" "reference FASTA in the index directory"

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
  --index "$test_data_dir/index" \
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
  --index "$test_data_dir/index" \
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
  --index "$test_data_dir/index" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/advanced.sam" \
  ---cpus 2 \
  --min_seed_length 15

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/advanced.sam" "advanced SAM output"
check_file_not_empty "$meta_temp_dir/advanced.sam" "advanced SAM output"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Index prefix unrelated to the reference file name ---
log "Starting TEST 4: BWA MEM with a custom index prefix"

# The prefix must be discovered from the directory contents, not derived from the
# name of the directory or of the reference FASTA.
log "Building BWA index with a custom prefix..."
mkdir -p "$test_data_dir/custom_index"
bwa-mem2 index -p "$test_data_dir/custom_index/my_genome_v2" "$test_data_dir/reference.fasta" >/dev/null 2>&1
check_file_exists "$test_data_dir/custom_index/my_genome_v2.bwt.2bit.64" "custom-prefixed index file"

log "Executing $meta_name with a custom index prefix..."
"$meta_executable" \
  --index "$test_data_dir/custom_index" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/custom_prefix.sam"

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/custom_prefix.sam" "custom-prefix SAM output"
check_file_not_empty "$meta_temp_dir/custom_prefix.sam" "custom-prefix SAM output"
check_file_contains "$meta_temp_dir/custom_prefix.sam" "@SQ" "custom-prefix SAM output"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Index built alongside the reference FASTA ---
log "Starting TEST 5: BWA MEM with an index built next to the reference FASTA"

# Without `-p` the index files are named after the FASTA, so the base name is
# `reference.fasta`. The prefix has to come from the index files themselves
# rather than from stripping a FASTA extension, and the reference FASTA sitting
# alongside the index must not confuse that.
log "Building BWA index next to the reference FASTA..."
mkdir -p "$test_data_dir/fasta_index"
cp "$test_data_dir/reference.fasta" "$test_data_dir/fasta_index/"
bwa-mem2 index "$test_data_dir/fasta_index/reference.fasta" >/dev/null 2>&1
check_file_exists "$test_data_dir/fasta_index/reference.fasta.bwt.2bit.64" "FASTA-named index file"

log "Executing $meta_name with an index built next to the reference FASTA..."
"$meta_executable" \
  --index "$test_data_dir/fasta_index" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/fasta_index.sam"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/fasta_index.sam" "FASTA-named-index SAM output"
check_file_not_empty "$meta_temp_dir/fasta_index.sam" "FASTA-named-index SAM output"
check_file_contains "$meta_temp_dir/fasta_index.sam" "@SQ" "FASTA-named-index SAM output"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Index directory without an index ---
log "Starting TEST 6: BWA MEM with an index directory that holds no index"

mkdir -p "$test_data_dir/empty_index"
cp "$test_data_dir/reference.fasta" "$test_data_dir/empty_index/"

log "Executing $meta_name with an index directory that holds no index..."
if "$meta_executable" \
  --index "$test_data_dir/empty_index" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/should_not_exist.sam" 2> "$meta_temp_dir/no_index.err"; then
  log_error "Expected a non-zero exit code for an index directory without an index"
  exit 1
fi

log "Validating TEST 6 outputs..."
check_file_contains "$meta_temp_dir/no_index.err" "No BWA-MEM2 index" "error message"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: Index directory holding more than one index ---
log "Starting TEST 7: BWA MEM with an ambiguous index directory"

mkdir -p "$test_data_dir/ambiguous_index"
cp "$test_data_dir/index"/reference.* "$test_data_dir/ambiguous_index/"
for ext in 0123 amb ann bwt.2bit.64 pac; do
  cp "$test_data_dir/index/reference.$ext" "$test_data_dir/ambiguous_index/second.$ext"
done

log "Executing $meta_name with an ambiguous index directory..."
if "$meta_executable" \
  --index "$test_data_dir/ambiguous_index" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/should_not_exist2.sam" 2> "$meta_temp_dir/ambiguous.err"; then
  log_error "Expected a non-zero exit code for an ambiguous index directory"
  exit 1
fi

log "Validating TEST 7 outputs..."
check_file_contains "$meta_temp_dir/ambiguous.err" "Multiple BWA-MEM2 indices" "error message"

log "✅ TEST 7 completed successfully"

# --- Test Case 8: A file passed where a directory is expected ---
log "Starting TEST 8: BWA MEM with a file instead of an index directory"

log "Executing $meta_name with a file instead of an index directory..."
if "$meta_executable" \
  --index "$test_data_dir/index/reference.bwt.2bit.64" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/should_not_exist3.sam" 2> "$meta_temp_dir/not_a_dir.err"; then
  log_error "Expected a non-zero exit code when --index is not a directory"
  exit 1
fi

log "Validating TEST 8 outputs..."
check_file_contains "$meta_temp_dir/not_a_dir.err" "must be a directory" "error message"

log "✅ TEST 8 completed successfully"

# --- Test Case 9: A broken index symlink alongside a valid index ---
log "Starting TEST 9: BWA MEM with a dangling index symlink in the index directory"

# A dangling '.bwt.2bit.64' symlink must not be counted as a second index: doing so
# would abort as ambiguous instead of using the index that is actually present.
mkdir -p "$test_data_dir/dangling_index"
cp "$test_data_dir/index"/reference.* "$test_data_dir/dangling_index/"
ln -s "does_not_exist.bwt.2bit.64" "$test_data_dir/dangling_index/dangling.bwt.2bit.64"

log "Executing $meta_name with a dangling index symlink present..."
"$meta_executable" \
  --index "$test_data_dir/dangling_index" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/dangling.sam"

log "Validating TEST 9 outputs..."
check_file_exists "$meta_temp_dir/dangling.sam" "dangling-symlink SAM output"
check_file_not_empty "$meta_temp_dir/dangling.sam" "dangling-symlink SAM output"
check_file_contains "$meta_temp_dir/dangling.sam" "@SQ" "dangling-symlink SAM output"

log "✅ TEST 9 completed successfully"

print_test_summary "All tests completed successfully"
