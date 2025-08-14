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
echo ">>> Test 1: Basic Bowtie2 indexing"
create_test_fasta "$TEMP_DIR/test_ref.fasta"

echo ">> Running bowtie2_build..."
"$meta_executable" \
  --input "$TEMP_DIR/test_ref.fasta" \
  --output "$TEMP_DIR/bt2_index"

echo ">> Checking output directory exists..."
assert_dir_exists "$TEMP_DIR/bt2_index"

echo ">> Checking all index files exist..."
assert_file_exists "$TEMP_DIR/bt2_index/test_ref.1.bt2"
assert_file_exists "$TEMP_DIR/bt2_index/test_ref.2.bt2"
assert_file_exists "$TEMP_DIR/bt2_index/test_ref.3.bt2"
assert_file_exists "$TEMP_DIR/bt2_index/test_ref.4.bt2"
assert_file_exists "$TEMP_DIR/bt2_index/test_ref.rev.1.bt2"
assert_file_exists "$TEMP_DIR/bt2_index/test_ref.rev.2.bt2"

echo ">> Checking index files are not empty..."
assert_file_not_empty "$TEMP_DIR/bt2_index/test_ref.1.bt2"
assert_file_not_empty "$TEMP_DIR/bt2_index/test_ref.2.bt2"
assert_file_not_empty "$TEMP_DIR/bt2_index/test_ref.3.bt2"
assert_file_not_empty "$TEMP_DIR/bt2_index/test_ref.4.bt2"
assert_file_not_empty "$TEMP_DIR/bt2_index/test_ref.rev.1.bt2"
assert_file_not_empty "$TEMP_DIR/bt2_index/test_ref.rev.2.bt2"

echo ">> OK: Basic indexing test passed."

# --- Test Case 2: Indexing with custom index name ---
echo ">>> Test 2: Bowtie2 indexing with custom index name"
create_test_fasta "$TEMP_DIR/genome.fasta"

"$meta_executable" \
  --input "$TEMP_DIR/genome.fasta" \
  --index_name "custom_genome" \
  --output "$TEMP_DIR/custom_index"

echo ">> Checking output directory exists..."
assert_dir_exists "$TEMP_DIR/custom_index"

echo ">> Checking all index files exist with custom name..."
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.1.bt2"
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.2.bt2"
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.3.bt2"
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.4.bt2"
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.rev.1.bt2"
assert_file_exists "$TEMP_DIR/custom_index/custom_genome.rev.2.bt2"

echo ">> Checking custom index files are not empty..."
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.1.bt2"
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.2.bt2"
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.3.bt2"
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.4.bt2"
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.rev.1.bt2"
assert_file_not_empty "$TEMP_DIR/custom_index/custom_genome.rev.2.bt2"

echo ">> OK: Custom index name test passed."

# --- Test Case 3: Indexing with --noref option ---
echo ">>> Test 3: Bowtie2 indexing with --noref"
create_test_fasta "$TEMP_DIR/ref_noref.fasta"

"$meta_executable" \
  --input "$TEMP_DIR/ref_noref.fasta" \
  --output "$TEMP_DIR/index_noref" \
  --noref

echo ">> Checking output directory exists..."
assert_dir_exists "$TEMP_DIR/index_noref"

echo ">> Checking forward index files exist..."
assert_file_exists "$TEMP_DIR/index_noref/ref_noref.1.bt2"
assert_file_exists "$TEMP_DIR/index_noref/ref_noref.2.bt2"
assert_file_exists "$TEMP_DIR/index_noref/ref_noref.rev.1.bt2"
assert_file_exists "$TEMP_DIR/index_noref/ref_noref.rev.2.bt2"

echo ">> Checking that .3/.4 index files are NOT created..."
if [ -f "$TEMP_DIR/index_noref/ref_noref.3.bt2" ] || [ -f "$TEMP_DIR/index_noref/ref_noref.4.bt2" ]; then
  echo "ERROR: .3/.4 index files should not exist with --noref option"
  exit 1
fi

echo ">> OK: --noref option test passed."

# --- Test Case 4: Error handling ---
echo ">>> Test 4: Error handling"

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
