#!/bin/bash

# Test Helper Functions for Biobox Components
# 
# This file provides standardized helper functions for component testing.
# Source this file in your test scripts with:
# source "$meta_resources_dir/test_helpers.sh"
#
# Usage examples:
#   log "Starting test execution"
#   check_file_exists "$output" "result file"
#   check_file_not_exists "$bam_file" "BAM file (disabled by default)"
#   create_test_fasta "$temp_dir/input.fasta" 3 50
#

#############################################
# Logging Functions
#############################################

# Log messages with timestamps and consistent formatting
log() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') [TEST] $*"
}

# Log informational messages (alias for log)
log_info() {
  log "$*"
}

# Log warning messages
log_warn() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') [WARN] $*"
}

# Log error messages
log_error() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] $*" >&2
}

#############################################
# File and Directory Validation Functions
#############################################

# Check if a file exists with descriptive logging
# Usage: check_file_exists "/path/to/file" "optional description"
check_file_exists() {
  local file_path="$1"
  local description="${2:-File}"
  
  if [[ -f "$file_path" ]]; then
    log "✓ Found $description: $file_path"
    return 0
  else
    log_error "✗ $description does not exist: $file_path"
    exit 1
  fi
}

# Check if a directory exists with descriptive logging
# Usage: check_dir_exists "/path/to/dir" "optional description"
check_dir_exists() {
  local dir_path="$1"
  local description="${2:-Directory}"
  
  if [[ -d "$dir_path" ]]; then
    log "✓ Found $description: $dir_path"
    return 0
  else
    log_error "✗ $description does not exist: $dir_path"
    exit 1
  fi
}

# Check if a file does NOT exist (useful for testing disabled features)
# Usage: check_file_not_exists "/path/to/file" "optional description"
check_file_not_exists() {
  local file_path="$1"
  local description="${2:-File}"
  
  if [[ ! -f "$file_path" ]]; then
    log "✓ Confirmed $description does not exist (as expected): $file_path"
    return 0
  else
    log_error "✗ $description exists but shouldn't: $file_path"
    exit 1
  fi
}

# Check if a directory does NOT exist (useful for testing disabled features)
# Usage: check_dir_not_exists "/path/to/dir" "optional description"
check_dir_not_exists() {
  local dir_path="$1"
  local description="${2:-Directory}"
  
  if [[ ! -d "$dir_path" ]]; then
    log "✓ Confirmed $description does not exist (as expected): $dir_path"
    return 0
  else
    log_error "✗ $description exists but shouldn't: $dir_path"
    exit 1
  fi
}

# Check if a file is not empty
# Usage: check_file_not_empty "/path/to/file" "optional description"
check_file_not_empty() {
  local file_path="$1"
  local description="${2:-File}"
  
  if [[ -s "$file_path" ]]; then
    log "✓ $description is not empty: $file_path"
    return 0
  else
    log_error "✗ $description is empty but shouldn't be: $file_path"
    exit 1
  fi
}

# Check if a file is empty
# Usage: check_file_empty "/path/to/file" "optional description"
check_file_empty() {
  local file_path="$1"
  local description="${2:-File}"
  
  if [[ ! -s "$file_path" ]]; then
    log "✓ $description is empty (as expected): $file_path"
    return 0
  else
    log_error "✗ $description is not empty but should be: $file_path"
    exit 1
  fi
}

#############################################
# Content Validation Functions
#############################################

# Check if a file contains specific text
# Usage: check_file_contains "/path/to/file" "search_text" "optional description"
check_file_contains() {
  local file_path="$1"
  local search_text="$2"
  local description="${3:-File}"
  
  if grep -q "$search_text" "$file_path" 2>/dev/null; then
    log "✓ $description contains expected text '$search_text': $file_path"
    return 0
  else
    log_error "✗ $description does not contain '$search_text': $file_path"
    exit 1
  fi
}

# Check if a file does NOT contain specific text
# Usage: check_file_not_contains "/path/to/file" "search_text" "optional description"
check_file_not_contains() {
  local file_path="$1"
  local search_text="$2"
  local description="${3:-File}"
  
  if ! grep -q "$search_text" "$file_path" 2>/dev/null; then
    log "✓ $description does not contain '$search_text' (as expected): $file_path"
    return 0
  else
    log_error "✗ $description contains '$search_text' but shouldn't: $file_path"
    exit 1
  fi
}

# Check if a file matches a regex pattern
# Usage: check_file_matches_regex "/path/to/file" "regex_pattern" "optional description"
check_file_matches_regex() {
  local file_path="$1"
  local regex_pattern="$2"
  local description="${3:-File}"
  
  if grep -qE "$regex_pattern" "$file_path" 2>/dev/null; then
    log "✓ $description matches expected pattern '$regex_pattern': $file_path"
    return 0
  else
    log_error "✗ $description does not match pattern '$regex_pattern': $file_path"
    exit 1
  fi
}

# Check if a file has the expected number of lines
# Usage: check_file_line_count "/path/to/file" expected_count "optional description"
check_file_line_count() {
  local file_path="$1"
  local expected_count="$2"
  local description="${3:-File}"
  
  local actual_count=$(wc -l < "$file_path" 2>/dev/null || echo "0")
  
  if [[ "$actual_count" -eq "$expected_count" ]]; then
    log "✓ $description has expected line count ($expected_count): $file_path"
    return 0
  else
    log_error "✗ $description has $actual_count lines, expected $expected_count: $file_path"
    exit 1
  fi
}

#############################################
# Test Data Generation Functions
#############################################

# Create a test FASTA file with specified sequences
# Usage: create_test_fasta "/path/to/output.fasta" [num_sequences] [sequence_length]
create_test_fasta() {
  local file_path="$1"
  local num_seqs="${2:-2}"
  local seq_length="${3:-64}"
  
  log "Creating test FASTA file with $num_seqs sequences of length $seq_length: $file_path"
  
  > "$file_path"  # Create empty file
  
  for i in $(seq 1 "$num_seqs"); do
    echo ">seq$i" >> "$file_path"
    # Generate random DNA sequence
    head -c "$seq_length" /dev/zero | tr '\0' 'A' | sed 's/A/ATCG/g' | head -c "$seq_length" >> "$file_path"
    echo >> "$file_path"
  done
  
  log "✓ Created test FASTA file: $file_path"
}

# Create a test FASTQ file with specified reads
# Usage: create_test_fastq "/path/to/output.fastq" [num_reads] [read_length]
create_test_fastq() {
  local file_path="$1"
  local num_reads="${2:-4}"
  local read_length="${3:-35}"
  
  log "Creating test FASTQ file with $num_reads reads of length $read_length: $file_path"
  
  > "$file_path"  # Create empty file
  
  for i in $(seq 1 "$num_reads"); do
    echo "@read$i" >> "$file_path"
    # Generate random DNA sequence
    head -c "$read_length" /dev/zero | tr '\0' 'A' | sed 's/A/ATCG/g' | head -c "$read_length" >> "$file_path"
    echo "+" >> "$file_path"
    # Generate quality scores (all high quality)
    head -c "$read_length" /dev/zero | tr '\0' 'I' >> "$file_path"
    echo >> "$file_path"
  done
  
  log "✓ Created test FASTQ file: $file_path"
}

# Create a test GTF file with basic gene annotations
# Usage: create_test_gtf "/path/to/output.gtf" [num_genes]
create_test_gtf() {
  local file_path="$1"
  local num_genes="${2:-3}"
  
  log "Creating test GTF file with $num_genes genes: $file_path"
  
  > "$file_path"  # Create empty file
  
  for i in $(seq 1 "$num_genes"); do
    local start=$((1000 * i))
    local end=$((start + 999))
    local chr="chr$((i % 22 + 1))"
    
    echo -e "${chr}\ttest\tgene\t${start}\t${end}\t.\t+\t.\tgene_id \"gene$i\"; gene_name \"GENE$i\"" >> "$file_path"
    echo -e "${chr}\ttest\ttranscript\t${start}\t${end}\t.\t+\t.\tgene_id \"gene$i\"; transcript_id \"transcript${i}\"; gene_name \"GENE$i\"" >> "$file_path"
    echo -e "${chr}\ttest\texon\t${start}\t$((start + 499))\t.\t+\t.\tgene_id \"gene$i\"; transcript_id \"transcript${i}\"; exon_number \"1\"" >> "$file_path"
    echo -e "${chr}\ttest\texon\t$((start + 500))\t${end}\t.\t+\t.\tgene_id \"gene$i\"; transcript_id \"transcript${i}\"; exon_number \"2\"" >> "$file_path"
  done
  
  log "✓ Created test GTF file: $file_path"
}

# Create a test GFF file with basic feature annotations
# Usage: create_test_gff "/path/to/output.gff" [num_features]
create_test_gff() {
  local file_path="$1"
  local num_features="${2:-3}"
  
  log "Creating test GFF file with $num_features features: $file_path"
  
  echo "##gff-version 3" > "$file_path"
  
  for i in $(seq 1 "$num_features"); do
    local start=$((1000 * i))
    local end=$((start + 999))
    local chr="chr$((i % 22 + 1))"
    
    echo -e "${chr}\ttest\tgene\t${start}\t${end}\t.\t+\t.\tID=gene$i;Name=GENE$i" >> "$file_path"
  done
  
  log "✓ Created test GFF file: $file_path"
}

# Create a test BED file with genomic intervals
# Usage: create_test_bed "/path/to/output.bed" [num_intervals]
create_test_bed() {
  local file_path="$1"
  local num_intervals="${2:-3}"
  
  log "Creating test BED file with $num_intervals intervals: $file_path"
  
  > "$file_path"  # Create empty file
  
  for i in $(seq 1 "$num_intervals"); do
    local start=$((1000 * i))
    local end=$((start + 999))
    local chr="chr$((i % 22 + 1))"
    
    echo -e "${chr}\t${start}\t${end}\tregion$i\t0\t+" >> "$file_path"
  done
  
  log "✓ Created test BED file: $file_path"
}

# Create a simple test CSV file
# Usage: create_test_csv "/path/to/output.csv" [num_rows]
create_test_csv() {
  local file_path="$1"
  local num_rows="${2:-5}"
  
  log "Creating test CSV file with $num_rows rows: $file_path"
  
  echo "id,name,value,category" > "$file_path"
  
  for i in $(seq 1 "$num_rows"); do
    echo "row$i,name$i,$((i * 10)),category$((i % 3 + 1))" >> "$file_path"
  done
  
  log "✓ Created test CSV file: $file_path"
}

# Create a simple test TSV file
# Usage: create_test_tsv "/path/to/output.tsv" [num_rows]
create_test_tsv() {
  local file_path="$1"
  local num_rows="${2:-5}"
  
  log "Creating test TSV file with $num_rows rows: $file_path"
  
  echo -e "id\tname\tvalue\tcategory" > "$file_path"
  
  for i in $(seq 1 "$num_rows"); do
    echo -e "row$i\tname$i\t$((i * 10))\tcategory$((i % 3 + 1))" >> "$file_path"
  done
  
  log "✓ Created test TSV file: $file_path"
}

#############################################
# Utility Functions
#############################################

# Setup test environment with recommended settings
setup_test_env() {
  # Enable strict error handling
  set -euo pipefail
  
  # Set up consistent locale for reproducible results
  export LC_ALL=C
  
  log "Test environment initialized with strict error handling"
  log "Using temporary directory: ${meta_temp_dir:-$PWD}"
}

# Clean up function (optional, since viash handles temp directory cleanup)
cleanup_test_env() {
  log "Test cleanup completed"
}

# Print test summary
print_test_summary() {
  local test_name="${1:-Test}"
  log "🎉 $test_name completed successfully!"
}

#############################################
# Example Usage
#############################################

# Example function showing how to use the helpers
example_test_usage() {
  log "=== Example Test Usage ==="
  
  # Setup
  setup_test_env
  
  # Create test data
  create_test_fasta "$meta_temp_dir/input.fasta" 3 50
  
  # Validate test data
  check_file_exists "$meta_temp_dir/input.fasta" "input FASTA file"
  check_file_not_empty "$meta_temp_dir/input.fasta" "input FASTA file"
  check_file_line_count "$meta_temp_dir/input.fasta" 6  # 3 sequences = 6 lines
  
  # Example tool execution (commented out)
  # "$meta_executable" --input "$meta_temp_dir/input.fasta" --output "$meta_temp_dir/output"
  
  # Validate outputs (examples)
  # check_file_exists "$meta_temp_dir/output.txt" "result file"
  # check_file_contains "$meta_temp_dir/output.txt" "expected_pattern" "result file"
  
  print_test_summary "Example test"
}
