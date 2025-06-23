#!/bin/bash

set -e

TEMP_DIR="$meta_temp_dir"

# --- Helper function to create a FASTQ file ---
create_fastq() {
  file_path="$1"
  num_records="$2"
  
  rm -f "$file_path"
  for i in $(seq 1 "$num_records"); do
    echo "@READ.${i} description" >> "$file_path"
    echo "GATTACA" >> "$file_path"
    echo "+" >> "$file_path"
    echo "FFFFFFF" >> "$file_path"
  done
}

# --- Test Case 1: Paired-End Subsampling with --record_count ---
echo ">>> Test 1: Paired-end subsampling with --record_count"
create_fastq "$TEMP_DIR/r1.fastq" 100
create_fastq "$TEMP_DIR/r2.fastq" 100

"$meta_executable" \
  --input_1 "$TEMP_DIR/r1.fastq" \
  --input_2 "$TEMP_DIR/r2.fastq" \
  --record_count 15 \
  --seed 42 \
  --output_1 "$TEMP_DIR/sub1.r1.fastq" \
  --output_2 "$TEMP_DIR/sub1.r2.fastq"

echo ">> Checking output files..."
if [ ! -f "$TEMP_DIR/sub1.r1.fastq" ]; then
  echo "FAIL: Subsampled R1 file was not created." && exit 1
fi
if [ ! -f "$TEMP_DIR/sub1.r2.fastq" ]; then
  echo "FAIL: Subsampled R2 file was not created." && exit 1
fi

# Each FASTQ record is 4 lines. 15 records * 4 lines/record = 60 lines.
line_count_r1=$(wc -l < "$TEMP_DIR/sub1.r1.fastq")
line_count_r2=$(wc -l < "$TEMP_DIR/sub1.r2.fastq")

if [ "$line_count_r1" -ne 60 ]; then
  echo "FAIL: R1 output has incorrect number of lines. Expected 60, got $line_count_r1." && exit 1
fi
if [ "$line_count_r2" -ne 60 ]; then
  echo "FAIL: R2 output has incorrect number of lines. Expected 60, got $line_count_r2." && exit 1
fi
echo ">> OK: Paired-end test with --record_count passed."

# --- Test Case 2: Single-End Subsampling with --probability and Gzipped Output ---
echo ">>> Test 2: Single-end subsampling with --probability and gzipped output"
create_fastq "$TEMP_DIR/r1.fastq" 500

"$meta_executable" \
  --input_1 "$TEMP_DIR/r1.fastq" \
  --probability 0.1 \
  --seed 42 \
  --output_1 "$TEMP_DIR/sub2.r1.fastq.gz"

echo ">> Checking gzipped output file..."
if [ ! -f "$TEMP_DIR/sub2.r1.fastq.gz" ]; then
  echo "FAIL: Gzipped subsampled file was not created." && exit 1
fi

# With a fixed seed, the number of records should be deterministic.
# NOTE: For fq v0.12.0, seed 42 and p=0.1 on 500 records yields 53 records. 53 * 4 = 212 lines.
gzipped_lines=$(gunzip -c "$TEMP_DIR/sub2.r1.fastq.gz" | wc -l)
if [ "$gzipped_lines" -ne 212 ]; then
  echo "FAIL: Gzipped output has incorrect number of lines. Expected 212, got $gzipped_lines." && exit 1
fi
echo ">> OK: Single-end test with --probability passed."


# --- Test Case 3: Mutually Exclusive Argument Check ---
echo ">>> Test 3: Expecting failure when both --record_count and --probability are provided"
set +e # Disable exit on error to catch the failure
"$meta_executable" \
  --input_1 "$TEMP_DIR/r1.fastq" \
  --record_count 10 \
  --probability 0.1 \
  --output_1 "$TEMP_DIR/sub3.r1.fastq"
exit_code=$?
set -e # Re-enable exit on error

if [ $exit_code -eq 0 ]; then
  echo "FAIL: Script should have failed when providing both count and probability."
  exit 1
else
  echo ">> OK: Script correctly failed as expected."
fi

echo ""
echo ">>> All tests finished successfully"
exit 0
