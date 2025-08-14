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
count_lines() {
  wc -l < "$1"
}
#############################################

# --- Test Case 1: Basic generation ---
echo ">>> Test 1: Basic FASTQ generation"

"$meta_executable" \
  --r1_dst "$TEMP_DIR/test_R1.fastq" \
  --r2_dst "$TEMP_DIR/test_R2.fastq" \
  --record_count 100 \
  --read_length 50

echo ">> Checking basic output..."
assert_file_exists "$TEMP_DIR/test_R1.fastq"
assert_file_exists "$TEMP_DIR/test_R2.fastq"
assert_file_not_empty "$TEMP_DIR/test_R1.fastq"
assert_file_not_empty "$TEMP_DIR/test_R2.fastq"

# Check FASTQ format
assert_file_contains "$TEMP_DIR/test_R1.fastq" "^@"
assert_file_contains "$TEMP_DIR/test_R2.fastq" "^@"
assert_file_contains "$TEMP_DIR/test_R1.fastq" "^+"
assert_file_contains "$TEMP_DIR/test_R2.fastq" "^+"

# Check record count (4 lines per record, so 100 records = 400 lines)
r1_lines=$(count_lines "$TEMP_DIR/test_R1.fastq")
r2_lines=$(count_lines "$TEMP_DIR/test_R2.fastq")

if [ "$r1_lines" -ne 400 ]; then
  echo "ERROR: Expected 400 lines in R1, got $r1_lines"
  exit 1
fi

if [ "$r2_lines" -ne 400 ]; then
  echo "ERROR: Expected 400 lines in R2, got $r2_lines"
  exit 1
fi

echo ">> OK: Basic generation test passed."

# --- Test Case 2: Generation with seed (reproducibility) ---
echo ">>> Test 2: Generation with seed for reproducibility"

"$meta_executable" \
  --r1_dst "$TEMP_DIR/seed1_R1.fastq" \
  --r2_dst "$TEMP_DIR/seed1_R2.fastq" \
  --record_count 50 \
  --read_length 30 \
  --seed 42

"$meta_executable" \
  --r1_dst "$TEMP_DIR/seed2_R1.fastq" \
  --r2_dst "$TEMP_DIR/seed2_R2.fastq" \
  --record_count 50 \
  --read_length 30 \
  --seed 42

echo ">> Checking reproducibility..."
if ! diff "$TEMP_DIR/seed1_R1.fastq" "$TEMP_DIR/seed2_R1.fastq" >/dev/null; then
  echo "ERROR: R1 files with same seed should be identical"
  exit 1
fi

if ! diff "$TEMP_DIR/seed1_R2.fastq" "$TEMP_DIR/seed2_R2.fastq" >/dev/null; then
  echo "ERROR: R2 files with same seed should be identical"
  exit 1
fi

echo ">> OK: Reproducibility test passed."

# --- Test Case 3: Gzipped output ---
echo ">>> Test 3: Gzipped output generation"

"$meta_executable" \
  --r1_dst "$TEMP_DIR/gzipped_R1.fastq.gz" \
  --r2_dst "$TEMP_DIR/gzipped_R2.fastq.gz" \
  --record_count 25 \
  --read_length 75

echo ">> Checking gzipped output..."
assert_file_exists "$TEMP_DIR/gzipped_R1.fastq.gz"
assert_file_exists "$TEMP_DIR/gzipped_R2.fastq.gz"
assert_file_not_empty "$TEMP_DIR/gzipped_R1.fastq.gz"
assert_file_not_empty "$TEMP_DIR/gzipped_R2.fastq.gz"

# Check that files are actually gzipped (check magic bytes)
if ! [ "$(head -c 2 "$TEMP_DIR/gzipped_R1.fastq.gz" | od -An -tx1)" = " 1f 8b" ]; then
  echo "ERROR: R1 file should be gzipped"
  exit 1
fi

if ! [ "$(head -c 2 "$TEMP_DIR/gzipped_R2.fastq.gz" | od -An -tx1)" = " 1f 8b" ]; then
  echo "ERROR: R2 file should be gzipped"
  exit 1
fi

# Test that gzipped files can be decompressed and contain valid FASTQ
gunzip -c "$TEMP_DIR/gzipped_R1.fastq.gz" | head -1 | grep -q "^@" || {
  echo "ERROR: Decompressed R1 file doesn't start with FASTQ header"
  exit 1
}

gunzip -c "$TEMP_DIR/gzipped_R2.fastq.gz" | head -1 | grep -q "^@" || {
  echo "ERROR: Decompressed R2 file doesn't start with FASTQ header"
  exit 1
}

echo ">> OK: Gzipped output test passed."

# --- Test Case 4: Different read lengths ---
echo ">>> Test 4: Custom read length"

"$meta_executable" \
  --r1_dst "$TEMP_DIR/custom_R1.fastq" \
  --r2_dst "$TEMP_DIR/custom_R2.fastq" \
  --record_count 10 \
  --read_length 150

echo ">> Checking custom read length..."
# Extract first sequence line and check length
seq_line=$(sed -n '2p' "$TEMP_DIR/custom_R1.fastq")
seq_length=${#seq_line}

if [ "$seq_length" -ne 150 ]; then
  echo "ERROR: Expected sequence length 150, got $seq_length"
  exit 1
fi

echo ">> OK: Custom read length test passed."

# --- Test Case 5: Default parameters ---
echo ">>> Test 5: Default parameters"

"$meta_executable" \
  --r1_dst "$TEMP_DIR/default_R1.fastq" \
  --r2_dst "$TEMP_DIR/default_R2.fastq"

echo ">> Checking default parameters..."
assert_file_exists "$TEMP_DIR/default_R1.fastq"
assert_file_exists "$TEMP_DIR/default_R2.fastq"
assert_file_not_empty "$TEMP_DIR/default_R1.fastq"
assert_file_not_empty "$TEMP_DIR/default_R2.fastq"

# Check default record count (10000 records = 40000 lines)
default_lines=$(count_lines "$TEMP_DIR/default_R1.fastq")
if [ "$default_lines" -ne 40000 ]; then
  echo "ERROR: Expected 40000 lines with default record count, got $default_lines"
  exit 1
fi

# Check default read length (101 bp)
default_seq_line=$(sed -n '2p' "$TEMP_DIR/default_R1.fastq")
default_seq_length=${#default_seq_line}

if [ "$default_seq_length" -ne 101 ]; then
  echo "ERROR: Expected default sequence length 101, got $default_seq_length"
  exit 1
fi

echo ">> OK: Default parameters test passed."

echo ">>> All tests passed!"
