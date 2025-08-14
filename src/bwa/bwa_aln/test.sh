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

# --- Helper function to create test reference ---
create_test_reference() {
  file_path="$1"
  
  cat << 'EOF' > "$file_path"
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF
}

# --- Helper function to create test FASTQ with shorter reads for BWA aln ---
create_test_fastq() {
  file_path="$1"
  read_prefix="$2"
  
  cat << EOF > "$file_path"
@${read_prefix}_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@${read_prefix}_2
CGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
}

#############################################

echo ">>> Setting up test reference and index"
create_test_reference "$TEMP_DIR/reference.fasta"

echo ">> Creating BWA index..."
cd "$TEMP_DIR"
bwa index reference.fasta
cd -

# --- Test Case 1: Basic BWA aln alignment ---
echo ">>> Test 1: BWA aln alignment"
create_test_fastq "$TEMP_DIR/reads.fastq" "read_aln"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --reads "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/output.sai"

echo ">> Checking output..."
assert_file_exists "$TEMP_DIR/output.sai"
assert_file_not_empty "$TEMP_DIR/output.sai"

# Check that it's a binary SAI file (should contain null bytes)
if grep -q "^@" "$TEMP_DIR/output.sai" 2>/dev/null; then
  echo "ERROR: SAI file appears to be text, should be binary"
  exit 1
fi

echo ">> OK: BWA aln alignment test passed."

# --- Test Case 2: BWA aln with custom parameters ---
echo ">>> Test 2: BWA aln with custom parameters"
create_test_fastq "$TEMP_DIR/reads_custom.fastq" "read_custom"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --reads "$TEMP_DIR/reads_custom.fastq" \
  --output "$TEMP_DIR/output_custom.sai" \
  --max_diff "0.05" \
  --seed_length 25 \
  --max_seed_diff 1 \
  --mismatch_penalty 4

echo ">> Checking custom parameter output..."
assert_file_exists "$TEMP_DIR/output_custom.sai"
assert_file_not_empty "$TEMP_DIR/output_custom.sai"

echo ">> OK: Custom parameter alignment test passed."

# --- Test Case 3: Error handling ---
echo ">>> Test 3: Error handling"

# Test with non-existent index
if "$meta_executable" \
  --index "$TEMP_DIR/nonexistent.fasta" \
  --reads "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/output_error.sai" 2>/dev/null; then
  echo "ERROR: Should have failed with non-existent index"
  exit 1
else
  echo ">> OK: Properly handled non-existent index error."
fi

echo ">>> All tests passed!"
