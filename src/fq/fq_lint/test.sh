#!/bin/bash

set -e

TEMP_DIR="$meta_temp_dir"

# --- Helper function to create FASTQ files ---
create_fastq() {
  file_path="$1"
  header_prefix="$2"
  num_records="$3"
  mismatch_quality="$4" # 'true' or 'false'
  
  rm -f "$file_path"
  for i in $(seq 1 "$num_records"); do
    seq="AATTGGCC"
    qual="FFFFFFFF"
    if [[ "$mismatch_quality" == "true" && "$i" -eq 2 ]]; then
      qual="FFFF" # Mismatched length for the second record
    fi
    echo "@${header_prefix}.${i} description" >> "$file_path"
    echo "$seq" >> "$file_path"
    echo "+" >> "$file_path"
    echo "$qual" >> "$file_path"
  done
}

# --- Test Case 1: Valid Paired-End FASTQ files ---
echo ">>> Test 1: Running on valid paired-end FASTQ files. Expecting success."
create_fastq "$TEMP_DIR/valid_r1.fastq" "PAIR" 10 "false"
create_fastq "$TEMP_DIR/valid_r2.fastq" "PAIR" 10 "false"

"$meta_executable" \
  --input_1 "$TEMP_DIR/valid_r1.fastq" \
  --input_2 "$TEMP_DIR/valid_r2.fastq"
echo ">> OK: fq lint succeeded on valid paired-end files."


# --- Test Case 2: Valid Single-End FASTQ file ---
echo ">>> Test 2: Running on a valid single-end FASTQ file. Expecting success."
"$meta_executable" \
  --input_1 "$TEMP_DIR/valid_r1.fastq"
echo ">> OK: fq lint succeeded on a valid single-end file."


# --- Test Case 3: Invalid Paired-End FASTQ (mismatched headers) ---
echo ">>> Test 3: Running on paired-end files with mismatched headers. Expecting failure."
create_fastq "$TEMP_DIR/mismatch_r1.fastq" "PAIR_A" 10 "false"
create_fastq "$TEMP_DIR/mismatch_r2.fastq" "PAIR_B" 10 "false"

# Disable exit on error temporarily to catch the expected failure
set +e
"$meta_executable" \
  --input_1 "$TEMP_DIR/mismatch_r1.fastq" \
  --input_2 "$TEMP_DIR/mismatch_r2.fastq"
exit_code=$?
set -e

if [ $exit_code -eq 0 ]; then
  echo ">> FAIL: fq lint should have failed on mismatched headers but succeeded."
  exit 1
else
  echo ">> OK: fq lint correctly failed on mismatched headers (Exit code: $exit_code)."
fi


# --- Test Case 4: Invalid Single-End FASTQ (sequence/quality length mismatch) ---
echo ">>> Test 4: Running on a single-end file with seq/qual length mismatch. Expecting failure."
create_fastq "$TEMP_DIR/bad_qual.fastq" "BAD" 5 "true"

set +e
"$meta_executable" \
  --input_1 "$TEMP_DIR/bad_qual.fastq"
exit_code=$?
set -e

if [ $exit_code -eq 0 ]; then
  echo ">> FAIL: fq lint should have failed on bad quality scores but succeeded."
  exit 1
else
  echo ">> OK: fq lint correctly failed on bad quality scores (Exit code: $exit_code)."
fi

# --- Test Case 5: Using --disable-validator to ignore mismatched headers ---
echo ">>> Test 5: Running on mismatched paired-end files but disabling validator P001. Expecting success."
# The validator for mismatched read names is P001 in `fq`.
"$meta_executable" \
  --input_1 "$TEMP_DIR/mismatch_r1.fastq" \
  --input_2 "$TEMP_DIR/mismatch_r2.fastq" \
  --disable_validator "P001"
echo ">> OK: fq lint succeeded when header mismatch validator was disabled."


echo ""
echo ">>> All tests finished successfully"
exit 0
