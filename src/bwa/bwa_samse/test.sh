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
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
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

echo ">>> Setting up SAI file using BWA aln"
create_test_fastq "$TEMP_DIR/reads.fastq" "read_samse"

echo ">> Running BWA aln to create SAI file..."
cd "$TEMP_DIR"
bwa aln reference.fasta reads.fastq > reads.sai
cd -

# --- Test Case 1: Basic BWA samse ---
echo ">>> Test 1: BWA samse basic functionality"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --sai "$TEMP_DIR/reads.sai" \
  --reads "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/output.sam"

echo ">> Checking output..."
assert_file_exists "$TEMP_DIR/output.sam"
assert_file_not_empty "$TEMP_DIR/output.sam"
assert_file_contains "$TEMP_DIR/output.sam" "@SQ"
assert_file_contains "$TEMP_DIR/output.sam" "@PG"

# Count alignment records
alignment_lines=$(grep -v "^@" "$TEMP_DIR/output.sam" | wc -l)
echo "Found $alignment_lines alignment records in SAM output."

echo ">> OK: BWA samse basic test passed."

# --- Test Case 2: BWA samse with read group ---
echo ">>> Test 2: BWA samse with read group"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --sai "$TEMP_DIR/reads.sai" \
  --reads "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/output_rg.sam" \
  --read_group "@RG\\tID:test\\tSM:sample1"

echo ">> Checking read group output..."
assert_file_exists "$TEMP_DIR/output_rg.sam"
assert_file_not_empty "$TEMP_DIR/output_rg.sam"
assert_file_contains "$TEMP_DIR/output_rg.sam" "@RG"
assert_file_contains "$TEMP_DIR/output_rg.sam" "ID:test"
assert_file_contains "$TEMP_DIR/output_rg.sam" "SM:sample1"

echo ">> OK: Read group test passed."

# --- Test Case 3: BWA samse with max occurrences ---
echo ">>> Test 3: BWA samse with max occurrences"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --sai "$TEMP_DIR/reads.sai" \
  --reads "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/output_maxocc.sam" \
  --max_occ 5

echo ">> Checking max occurrences output..."
assert_file_exists "$TEMP_DIR/output_maxocc.sam"
assert_file_not_empty "$TEMP_DIR/output_maxocc.sam"

echo ">> OK: Max occurrences test passed."

# --- Test Case 4: Error handling ---
echo ">>> Test 4: Error handling"

# Test with non-existent SAI file
if "$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --sai "$TEMP_DIR/nonexistent.sai" \
  --reads "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/output_error.sam" 2>/dev/null; then
  echo "ERROR: Should have failed with non-existent SAI file"
  exit 1
else
  echo ">> OK: Properly handled non-existent SAI file error."
fi

echo ">>> All tests passed!"
