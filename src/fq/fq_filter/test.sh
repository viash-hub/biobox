#!/bin/bash

set -e

TEMP_DIR="$meta_temp_dir"

# --- Helper function to create a FASTQ file ---
create_fastq() {
  file_path="$1"
  
  rm -f "$file_path"
  cat << 'EOF' >> "$file_path"
@READ.1 description
GATTACA
+
FFFFFFF
@READ.2 description
ATGCATG
+
HHHHHHH
@READ.3 description
CCCGGGG
+
IIIIIII
@READ.4 description
GATTACA
+
JJJJJJJ
@READ.5 description
TTTAAAA
+
KKKKKKK
EOF
}

# --- Helper function to create names file ---
create_names_file() {
  file_path="$1"
  
  rm -f "$file_path"
  cat << 'EOF' >> "$file_path"
READ.1
READ.3
READ.5
EOF
}

# --- Test Case 1: Filter by names ---
echo ">>> Test 1: Filter FASTQ by record names"
create_fastq "$TEMP_DIR/input.fastq"
create_names_file "$TEMP_DIR/names.txt"

"$meta_executable" \
  --input "$TEMP_DIR/input.fastq" \
  --output "$TEMP_DIR/filtered_names.fastq" \
  --names "$TEMP_DIR/names.txt"

echo ">> Checking filtered output..."
if [ ! -f "$TEMP_DIR/filtered_names.fastq" ]; then
  echo "FAIL: Filtered file was not created." && exit 1
fi

# Should have 3 records (READ.1, READ.3, READ.5) = 12 lines
line_count=$(wc -l < "$TEMP_DIR/filtered_names.fastq")
if [ "$line_count" -ne 12 ]; then
  echo "FAIL: Filtered output has incorrect number of lines. Expected 12, got $line_count." && exit 1
fi

# Check that the correct reads are present
if ! grep -q "READ.1" "$TEMP_DIR/filtered_names.fastq"; then
  echo "FAIL: READ.1 not found in filtered output." && exit 1
fi
if ! grep -q "READ.3" "$TEMP_DIR/filtered_names.fastq"; then
  echo "FAIL: READ.3 not found in filtered output." && exit 1
fi
if ! grep -q "READ.5" "$TEMP_DIR/filtered_names.fastq"; then
  echo "FAIL: READ.5 not found in filtered output." && exit 1
fi
if grep -q "READ.2" "$TEMP_DIR/filtered_names.fastq"; then
  echo "FAIL: READ.2 should not be in filtered output." && exit 1
fi
if grep -q "READ.4" "$TEMP_DIR/filtered_names.fastq"; then
  echo "FAIL: READ.4 should not be in filtered output." && exit 1
fi

echo ">> OK: Names filtering test passed."

# --- Test Case 2: Filter by sequence pattern ---
echo ">>> Test 2: Filter FASTQ by sequence pattern"
create_fastq "$TEMP_DIR/input2.fastq"

"$meta_executable" \
  --input "$TEMP_DIR/input2.fastq" \
  --output "$TEMP_DIR/filtered_pattern.fastq" \
  --sequence_pattern "GATTACA"

echo ">> Checking pattern filtered output..."
if [ ! -f "$TEMP_DIR/filtered_pattern.fastq" ]; then
  echo "FAIL: Pattern filtered file was not created." && exit 1
fi

# Should have 2 records (READ.1 and READ.4 have GATTACA) = 8 lines
line_count=$(wc -l < "$TEMP_DIR/filtered_pattern.fastq")
if [ "$line_count" -ne 8 ]; then
  echo "FAIL: Pattern filtered output has incorrect number of lines. Expected 8, got $line_count." && exit 1
fi

# Check that the correct reads are present
if ! grep -q "READ.1" "$TEMP_DIR/filtered_pattern.fastq"; then
  echo "FAIL: READ.1 not found in pattern filtered output." && exit 1
fi
if ! grep -q "READ.4" "$TEMP_DIR/filtered_pattern.fastq"; then
  echo "FAIL: READ.4 not found in pattern filtered output." && exit 1
fi
if grep -q "READ.2" "$TEMP_DIR/filtered_pattern.fastq"; then
  echo "FAIL: READ.2 should not be in pattern filtered output." && exit 1
fi

echo ">> OK: Pattern filtering test passed."

# --- Test Case 3: Test with gzipped output ---
echo ">>> Test 3: Filter with gzipped output"
create_fastq "$TEMP_DIR/input3.fastq"

"$meta_executable" \
  --input "$TEMP_DIR/input3.fastq" \
  --output "$TEMP_DIR/filtered.fastq.gz" \
  --sequence_pattern "ATG"

echo ">> Checking gzipped output..."
if [ ! -f "$TEMP_DIR/filtered.fastq.gz" ]; then
  echo "FAIL: Gzipped filtered file was not created." && exit 1
fi

# Should have 1 record (READ.2 has ATGCATG) = 4 lines
gzipped_lines=$(gunzip -c "$TEMP_DIR/filtered.fastq.gz" | wc -l)
if [ "$gzipped_lines" -ne 4 ]; then
  echo "FAIL: Gzipped output has incorrect number of lines. Expected 4, got $gzipped_lines." && exit 1
fi

echo ">> OK: Gzipped output test passed."

# --- Test Case 4: Test error when no filtering options provided ---
echo ">>> Test 4: Expecting failure when no filtering options are provided"
set +e # Disable exit on error to catch the failure
"$meta_executable" \
  --input "$TEMP_DIR/input.fastq" \
  --output "$TEMP_DIR/error_test.fastq"
exit_code=$?
set -e # Re-enable exit on error

if [ $exit_code -eq 0 ]; then
  echo "FAIL: Script should have failed when no filtering options are provided."
  exit 1
else
  echo ">> OK: Script correctly failed as expected."
fi

echo ""
echo ">>> All tests finished successfully"
exit 0
