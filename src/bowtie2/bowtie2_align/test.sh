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

# --- Helper functions to create test data ---
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

create_test_fastq() {
  file_path="$1"
  
  cat << 'EOF' > "$file_path"
@read1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read4
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
}

create_test_fastq_paired() {
  file_path_r1="$1"
  file_path_r2="$2"
  
  cat << 'EOF' > "$file_path_r1"
@read1/1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2/1
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

  cat << 'EOF' > "$file_path_r2"
@read1/2
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2/2
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
}

# --- Prepare test data ---
echo ">>> Preparing test data"
create_test_fasta "$TEMP_DIR/test_ref.fasta"

# Build index using bowtie2-build
echo ">> Building Bowtie2 index..."
mkdir -p "$TEMP_DIR/index"
bowtie2-build "$TEMP_DIR/test_ref.fasta" "$TEMP_DIR/index/test_ref" > /dev/null 2>&1

create_test_fastq "$TEMP_DIR/reads.fastq"
create_test_fastq_paired "$TEMP_DIR/reads_R1.fastq" "$TEMP_DIR/reads_R2.fastq"

# --- Test Case 1: Single-end alignment ---
echo ">>> Test 1: Single-end alignment"

"$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --unpaired "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/single_end.sam"

echo ">> Checking single-end output..."
assert_file_exists "$TEMP_DIR/single_end.sam"
assert_file_not_empty "$TEMP_DIR/single_end.sam"

# Check for SAM header
if ! grep -q "^@HD" "$TEMP_DIR/single_end.sam"; then
  echo "ERROR: SAM header not found"
  exit 1
fi

echo ">> OK: Single-end alignment test passed."

# --- Test Case 2: Paired-end alignment ---
echo ">>> Test 2: Paired-end alignment"

"$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --mate1 "$TEMP_DIR/reads_R1.fastq" \
  --mate2 "$TEMP_DIR/reads_R2.fastq" \
  --output "$TEMP_DIR/paired_end.sam"

echo ">> Checking paired-end output..."
assert_file_exists "$TEMP_DIR/paired_end.sam"
assert_file_not_empty "$TEMP_DIR/paired_end.sam"

# Check for SAM header
if ! grep -q "^@HD" "$TEMP_DIR/paired_end.sam"; then
  echo "ERROR: SAM header not found"
  exit 1
fi

echo ">> OK: Paired-end alignment test passed."

# --- Test Case 3: Alignment with preset ---
echo ">>> Test 3: Alignment with very-fast preset"

"$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --unpaired "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/very_fast.sam" \
  --very_fast

echo ">> Checking very-fast preset output..."
assert_file_exists "$TEMP_DIR/very_fast.sam"
assert_file_not_empty "$TEMP_DIR/very_fast.sam"

echo ">> OK: Very-fast preset test passed."

# --- Test Case 4: Local alignment ---
echo ">>> Test 4: Local alignment"

"$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --unpaired "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/local.sam" \
  --local

echo ">> Checking local alignment output..."
assert_file_exists "$TEMP_DIR/local.sam"
assert_file_not_empty "$TEMP_DIR/local.sam"

echo ">> OK: Local alignment test passed."

# --- Test Case 5: Alignment with unaligned output ---
echo ">>> Test 5: Alignment with unaligned output"

"$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --unpaired "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/with_unaligned.sam" \
  --un "$TEMP_DIR/unaligned.fastq"

echo ">> Checking output with unaligned reads..."
assert_file_exists "$TEMP_DIR/with_unaligned.sam"
assert_file_not_empty "$TEMP_DIR/with_unaligned.sam"

# Unaligned file may or may not exist depending on whether all reads aligned
if [ -f "$TEMP_DIR/unaligned.fastq" ]; then
  echo ">> Unaligned reads file created"
else
  echo ">> No unaligned reads (all reads aligned successfully)"
fi

echo ">> OK: Unaligned output test passed."

# --- Test Case 6: Error handling ---
echo ">>> Test 6: Error handling"

# Test with non-existent index
if "$meta_executable" \
  --index "$TEMP_DIR/nonexistent_index" \
  --unpaired "$TEMP_DIR/reads.fastq" \
  --output "$TEMP_DIR/error.sam" 2>/dev/null; then
  echo "ERROR: Should have failed with non-existent index"
  exit 1
else
  echo ">> OK: Properly handled non-existent index error."
fi

# Test with missing input files
if "$meta_executable" \
  --index "$TEMP_DIR/index/test_ref" \
  --unpaired "$TEMP_DIR/nonexistent.fastq" \
  --output "$TEMP_DIR/error.sam" 2>/dev/null; then
  echo "ERROR: Should have failed with non-existent input file"
  exit 1
else
  echo ">> OK: Properly handled non-existent input file error."
fi

echo ">>> All tests passed!"
