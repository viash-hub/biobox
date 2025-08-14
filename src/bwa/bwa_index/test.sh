#!/bin/bash

set -e

TEMP_DIR="$meta_temp_dir"

#############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_dir_exists() {
  [ -d "$1" ] || { echo "Directory '$1' does not exist" && exit 1; }
}
#############################################

# --- Helper function to create test FASTA ---
create_test_fasta() {
  file_path="$1"
  
  cat << 'EOF' > "$file_path"
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
EOF
}

# --- Test Case 1: Basic indexing ---
echo ">>> Test 1: Basic BWA indexing"
create_test_fasta "$TEMP_DIR/test_ref.fasta"

echo ">> Running bwa_index..."
"$meta_executable" \
  --input "$TEMP_DIR/test_ref.fasta" \
  --output "$TEMP_DIR/bwa_index"

echo ">> Checking output directory exists..."
assert_dir_exists "$TEMP_DIR/bwa_index"

echo ">> Checking all index files exist..."
assert_file_exists "$TEMP_DIR/bwa_index/test_ref.amb"
assert_file_exists "$TEMP_DIR/bwa_index/test_ref.ann"
assert_file_exists "$TEMP_DIR/bwa_index/test_ref.bwt"
assert_file_exists "$TEMP_DIR/bwa_index/test_ref.pac"
assert_file_exists "$TEMP_DIR/bwa_index/test_ref.sa"

echo ">> Checking index files are not empty..."
assert_file_not_empty "$TEMP_DIR/bwa_index/test_ref.amb"
assert_file_not_empty "$TEMP_DIR/bwa_index/test_ref.ann"
assert_file_not_empty "$TEMP_DIR/bwa_index/test_ref.bwt"
assert_file_not_empty "$TEMP_DIR/bwa_index/test_ref.pac"
assert_file_not_empty "$TEMP_DIR/bwa_index/test_ref.sa"

echo ">> OK: Basic indexing test passed."

# --- Test Case 2: Indexing with custom prefix ---
echo ">>> Test 2: BWA indexing with custom prefix"
create_test_fasta "$TEMP_DIR/genome.fasta"

"$meta_executable" \
  --input "$TEMP_DIR/genome.fasta" \
  --prefix "custom_genome" \
  --output "$TEMP_DIR/custom_index"

echo ">> Checking output directory exists..."
assert_dir_exists "$TEMP_DIR/custom_index"

echo ">> Checking all index files exist with custom prefix..."
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.amb"
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.ann"
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.bwt"
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.pac"
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.sa"

echo ">> Checking custom prefix index files are not empty..."
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.amb"
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.ann"
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.bwt"
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.pac"
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.sa"

echo ">> OK: Custom prefix indexing test passed."

# --- Test Case 3: Error handling ---
echo ">>> Test 3: Error handling"

# Test with non-existent input file
if "$meta_executable" \
  --input "$TEMP_DIR/nonexistent.fasta" \
  --output "$TEMP_DIR/error_index" 2>/dev/null; then
  echo "ERROR: Should have failed with non-existent input file"
  exit 1
else
  echo ">> OK: Properly handled non-existent input file error."
fi

echo ">>> All tests passed!"
