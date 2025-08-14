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
  read_suffix="$3"
  
  cat << EOF > "$file_path"
@${read_prefix}_1_${read_suffix}
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@${read_prefix}_2_${read_suffix}
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

echo ">>> Setting up SAI files using BWA aln"
create_test_fastq "$TEMP_DIR/reads_R1.fastq" "read_sampe" "R1"
create_test_fastq "$TEMP_DIR/reads_R2.fastq" "read_sampe" "R2"

echo ">> Running BWA aln to create SAI files..."
cd "$TEMP_DIR"
bwa aln reference.fasta reads_R1.fastq > reads_R1.sai
bwa aln reference.fasta reads_R2.fastq > reads_R2.sai
cd -

# --- Test Case 1: Basic BWA sampe ---
echo ">>> Test 1: BWA sampe basic functionality"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --sai1 "$TEMP_DIR/reads_R1.sai" \
  --sai2 "$TEMP_DIR/reads_R2.sai" \
  --reads1 "$TEMP_DIR/reads_R1.fastq" \
  --reads2 "$TEMP_DIR/reads_R2.fastq" \
  --output "$TEMP_DIR/output.sam"

echo ">> Checking output..."
assert_file_exists "$TEMP_DIR/output.sam"
assert_file_not_empty "$TEMP_DIR/output.sam"
assert_file_contains "$TEMP_DIR/output.sam" "@SQ"
assert_file_contains "$TEMP_DIR/output.sam" "@PG"

# Count alignment records
alignment_lines=$(grep -v "^@" "$TEMP_DIR/output.sam" | wc -l)
echo "Found $alignment_lines alignment records in SAM output."

echo ">> OK: BWA sampe basic test passed."

# --- Test Case 2: BWA sampe with read group ---
echo ">>> Test 2: BWA sampe with read group"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --sai1 "$TEMP_DIR/reads_R1.sai" \
  --sai2 "$TEMP_DIR/reads_R2.sai" \
  --reads1 "$TEMP_DIR/reads_R1.fastq" \
  --reads2 "$TEMP_DIR/reads_R2.fastq" \
  --output "$TEMP_DIR/output_rg.sam" \
  --read_group "@RG\\tID:test\\tSM:sample1"

echo ">> Checking read group output..."
assert_file_exists "$TEMP_DIR/output_rg.sam"
assert_file_not_empty "$TEMP_DIR/output_rg.sam"
assert_file_contains "$TEMP_DIR/output_rg.sam" "@RG"
assert_file_contains "$TEMP_DIR/output_rg.sam" "ID:test"
assert_file_contains "$TEMP_DIR/output_rg.sam" "SM:sample1"

echo ">> OK: Read group test passed."

# --- Test Case 3: BWA sampe with custom parameters ---
echo ">>> Test 3: BWA sampe with custom parameters"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --sai1 "$TEMP_DIR/reads_R1.sai" \
  --sai2 "$TEMP_DIR/reads_R2.sai" \
  --reads1 "$TEMP_DIR/reads_R1.fastq" \
  --reads2 "$TEMP_DIR/reads_R2.fastq" \
  --output "$TEMP_DIR/output_custom.sam" \
  --max_insert_size 300 \
  --max_hits_paired 5 \
  --max_hits_discordant 15

echo ">> Checking custom parameters output..."
assert_file_exists "$TEMP_DIR/output_custom.sam"
assert_file_not_empty "$TEMP_DIR/output_custom.sam"

echo ">> OK: Custom parameters test passed."

# --- Test Case 4: BWA sampe with algorithm options ---
echo ">>> Test 4: BWA sampe with algorithm options"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --sai1 "$TEMP_DIR/reads_R1.sai" \
  --sai2 "$TEMP_DIR/reads_R2.sai" \
  --reads1 "$TEMP_DIR/reads_R1.fastq" \
  --reads2 "$TEMP_DIR/reads_R2.fastq" \
  --output "$TEMP_DIR/output_algo.sam" \
  --disable_smith_waterman \
  --disable_insert_size_estimate

echo ">> Checking algorithm options output..."
assert_file_exists "$TEMP_DIR/output_algo.sam"
assert_file_not_empty "$TEMP_DIR/output_algo.sam"

echo ">> OK: Algorithm options test passed."

# --- Test Case 5: Error handling ---
echo ">>> Test 5: Error handling"

# Test with non-existent SAI file
if "$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --sai1 "$TEMP_DIR/nonexistent.sai" \
  --sai2 "$TEMP_DIR/reads_R2.sai" \
  --reads1 "$TEMP_DIR/reads_R1.fastq" \
  --reads2 "$TEMP_DIR/reads_R2.fastq" \
  --output "$TEMP_DIR/output_error.sam" 2>/dev/null; then
  echo "ERROR: Should have failed with non-existent SAI file"
  exit 1
else
  echo ">> OK: Properly handled non-existent SAI file error."
fi

echo ">>> All tests passed!"
