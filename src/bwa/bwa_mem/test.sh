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

# --- Helper function to create test FASTQ ---
create_test_fastq() {
  file_path="$1"
  read_prefix="$2"
  
  cat << EOF > "$file_path"
@${read_prefix}_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@${read_prefix}_2
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
}

#############################################

echo ">>> Setting up test reference and index"
create_test_reference "$TEMP_DIR/reference.fasta"

echo ">> Creating BWA index..."
cd "$TEMP_DIR"
bwa index reference.fasta
cd -

# --- Test Case 1: Single-end alignment ---
echo ">>> Test 1: Single-end BWA-MEM alignment"
create_test_fastq "$TEMP_DIR/reads_SE.fastq" "read_SE"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --reads1 "$TEMP_DIR/reads_SE.fastq" \
  --output "$TEMP_DIR/output_SE.sam"

echo ">> Checking single-end output..."
assert_file_exists "$TEMP_DIR/output_SE.sam"
assert_file_not_empty "$TEMP_DIR/output_SE.sam"
assert_file_contains "$TEMP_DIR/output_SE.sam" "@SQ"
assert_file_contains "$TEMP_DIR/output_SE.sam" "@PG"
# Count alignment records
alignment_lines=$(grep -vc "^@" "$TEMP_DIR/output_SE.sam")
echo "Found $alignment_lines alignment records in single-end output."

echo ">> OK: Single-end alignment test passed."

# --- Test Case 2: Paired-end alignment ---
echo ">>> Test 2: Paired-end BWA-MEM alignment"
create_test_fastq "$TEMP_DIR/reads_R1.fastq" "read_PE"
create_test_fastq "$TEMP_DIR/reads_R2.fastq" "read_PE"

"$meta_executable" \
  --index "$TEMP_DIR/reference.fasta" \
  --reads1 "$TEMP_DIR/reads_R1.fastq" \
  --reads2 "$TEMP_DIR/reads_R2.fastq" \
  --output "$TEMP_DIR/output_PE.sam"

echo ">> Checking paired-end output..."
assert_file_exists "$TEMP_DIR/output_PE.sam"
assert_file_not_empty "$TEMP_DIR/output_PE.sam"
assert_file_contains "$TEMP_DIR/output_PE.sam" "@SQ"
assert_file_contains "$TEMP_DIR/output_PE.sam" "@PG"
# Count alignment records
alignment_lines=$(grep -vc "^@" "$TEMP_DIR/output_PE.sam")
echo "Found $alignment_lines alignment records in paired-end output."

echo ">> OK: Paired-end alignment test passed."

echo ">>> All tests passed!"
