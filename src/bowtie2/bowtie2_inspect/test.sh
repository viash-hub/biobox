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

# --- Prepare test data ---
echo ">>> Preparing test data"
create_test_fasta "$TEMP_DIR/test_ref.fasta"

# Build index using bowtie2-build
echo ">> Building Bowtie2 index..."
mkdir -p "$TEMP_DIR/index"
bowtie2-build "$TEMP_DIR/test_ref.fasta" "$TEMP_DIR/index/test_ref" > /dev/null 2>&1

# --- Test Case 1: Default FASTA output ---
echo ">>> Test 1: Default FASTA output"

"$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --output "$TEMP_DIR/sequences.fasta"

echo ">> Checking FASTA output..."
assert_file_exists "$TEMP_DIR/sequences.fasta"
assert_file_not_empty "$TEMP_DIR/sequences.fasta"

# Check for FASTA format
if ! grep -q "^>" "$TEMP_DIR/sequences.fasta"; then
  echo "ERROR: FASTA headers not found"
  exit 1
fi

if ! grep -q "chr1" "$TEMP_DIR/sequences.fasta"; then
  echo "ERROR: Expected sequence chr1 not found"
  exit 1
fi

if ! grep -q "chr2" "$TEMP_DIR/sequences.fasta"; then
  echo "ERROR: Expected sequence chr2 not found"
  exit 1
fi

echo ">> OK: Default FASTA output test passed."

# --- Test Case 2: Summary output ---
echo ">>> Test 2: Summary output"

"$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --summary \
  --output "$TEMP_DIR/summary.txt"

echo ">> Checking summary output..."
assert_file_exists "$TEMP_DIR/summary.txt"
assert_file_not_empty "$TEMP_DIR/summary.txt"

# Check for summary content
if ! grep -q "Sequence" "$TEMP_DIR/summary.txt"; then
  echo "ERROR: Summary should contain 'Sequence' information"
  exit 1
fi

echo ">> OK: Summary output test passed."

# --- Test Case 3: Names only output ---
echo ">>> Test 3: Names only output"

"$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --names \
  --output "$TEMP_DIR/names.txt"

echo ">> Checking names output..."
assert_file_exists "$TEMP_DIR/names.txt"
assert_file_not_empty "$TEMP_DIR/names.txt"

# Check for sequence names
if ! grep -q "chr1" "$TEMP_DIR/names.txt"; then
  echo "ERROR: Expected sequence name chr1 not found"
  exit 1
fi

if ! grep -q "chr2" "$TEMP_DIR/names.txt"; then
  echo "ERROR: Expected sequence name chr2 not found"
  exit 1
fi

# Make sure it doesn't contain actual sequences (no ATCG pattern)
if grep -q "ATCG" "$TEMP_DIR/names.txt"; then
  echo "ERROR: Names output should not contain sequence data"
  exit 1
fi

echo ">> OK: Names only output test passed."

# --- Test Case 4: Custom line width ---
echo ">>> Test 4: Custom line width"

"$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --across 20 \
  --output "$TEMP_DIR/wide_sequences.fasta"

echo ">> Checking custom line width output..."
assert_file_exists "$TEMP_DIR/wide_sequences.fasta"
assert_file_not_empty "$TEMP_DIR/wide_sequences.fasta"

# Check for FASTA format
if ! grep -q "^>" "$TEMP_DIR/wide_sequences.fasta"; then
  echo "ERROR: FASTA headers not found"
  exit 1
fi

echo ">> OK: Custom line width test passed."

# --- Test Case 5: Error handling ---
echo ">>> Test 5: Error handling"

# Test with non-existent index
if "$meta_executable" \
  --index "$TEMP_DIR/nonexistent_index" \
  --output "$TEMP_DIR/error.txt" 2>/dev/null; then
  echo "ERROR: Should have failed with non-existent index"
  exit 1
else
  echo ">> OK: Properly handled non-existent index error."
fi

echo ">>> All tests passed!"
