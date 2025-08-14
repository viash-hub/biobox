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
#############################################

# --- Helper function to create a test FASTA file ---
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
  --output_amb "$TEMP_DIR/test_ref.fasta.amb" \
  --output_ann "$TEMP_DIR/test_ref.fasta.ann" \
  --output_bwt "$TEMP_DIR/test_ref.fasta.bwt" \
  --output_pac "$TEMP_DIR/test_ref.fasta.pac" \
  --output_sa "$TEMP_DIR/test_ref.fasta.sa"

echo ">> Checking all index files exist..."
assert_file_exists "$TEMP_DIR/test_ref.fasta.amb"
assert_file_exists "$TEMP_DIR/test_ref.fasta.ann"
assert_file_exists "$TEMP_DIR/test_ref.fasta.bwt"
assert_file_exists "$TEMP_DIR/test_ref.fasta.pac"
assert_file_exists "$TEMP_DIR/test_ref.fasta.sa"

echo ">> Checking index files are not empty..."
assert_file_not_empty "$TEMP_DIR/test_ref.fasta.amb"
assert_file_not_empty "$TEMP_DIR/test_ref.fasta.ann"
assert_file_not_empty "$TEMP_DIR/test_ref.fasta.bwt"
assert_file_not_empty "$TEMP_DIR/test_ref.fasta.pac"
assert_file_not_empty "$TEMP_DIR/test_ref.fasta.sa"

echo ">> OK: Basic indexing test passed."

# --- Test Case 2: Indexing with custom prefix ---
echo ">>> Test 2: BWA indexing with custom prefix"
create_test_fasta "$TEMP_DIR/genome.fasta"

"$meta_executable" \
  --input "$TEMP_DIR/genome.fasta" \
  --prefix "$TEMP_DIR/custom_index" \
  --output_amb "$TEMP_DIR/custom_index.amb" \
  --output_ann "$TEMP_DIR/custom_index.ann" \
  --output_bwt "$TEMP_DIR/custom_index.bwt" \
  --output_pac "$TEMP_DIR/custom_index.pac" \
  --output_sa "$TEMP_DIR/custom_index.sa"

echo ">> Checking all index files exist with custom prefix..."
assert_file_exists "$TEMP_DIR/custom_index.amb"
assert_file_exists "$TEMP_DIR/custom_index.ann"
assert_file_exists "$TEMP_DIR/custom_index.bwt"
assert_file_exists "$TEMP_DIR/custom_index.pac"
assert_file_exists "$TEMP_DIR/custom_index.sa"

echo ">> Checking index files are not empty..."
assert_file_not_empty "$TEMP_DIR/custom_index.amb"
assert_file_not_empty "$TEMP_DIR/custom_index.ann"
assert_file_not_empty "$TEMP_DIR/custom_index.bwt"
assert_file_not_empty "$TEMP_DIR/custom_index.pac"
assert_file_not_empty "$TEMP_DIR/custom_index.sa"

echo ">> OK: Custom prefix test passed."

# --- Test Case 3: Test with 64-bit naming ---
echo ">>> Test 3: BWA indexing with 64-bit naming"
create_test_fasta "$TEMP_DIR/ref64.fasta"

"$meta_executable" \
  --input "$TEMP_DIR/ref64.fasta" \
  --use_64bit_names \
  --output_amb "$TEMP_DIR/ref64.fasta.64.amb" \
  --output_ann "$TEMP_DIR/ref64.fasta.64.ann" \
  --output_bwt "$TEMP_DIR/ref64.fasta.64.bwt" \
  --output_pac "$TEMP_DIR/ref64.fasta.64.pac" \
  --output_sa "$TEMP_DIR/ref64.fasta.64.sa"

echo ">> Checking all 64-bit index files exist..."
assert_file_exists "$TEMP_DIR/ref64.fasta.64.amb"
assert_file_exists "$TEMP_DIR/ref64.fasta.64.ann"
assert_file_exists "$TEMP_DIR/ref64.fasta.64.bwt"
assert_file_exists "$TEMP_DIR/ref64.fasta.64.pac"
assert_file_exists "$TEMP_DIR/ref64.fasta.64.sa"

echo ">> OK: 64-bit naming test passed."

echo ""
echo ">>> All tests finished successfully"
exit 0
