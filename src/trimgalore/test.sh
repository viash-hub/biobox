#!/bin/bash

set -e
set -eo pipefail

# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_doesnt_exist() {
  [ ! -f "$1" ] || { echo "File '$1' exists but shouldn't" && exit 1; }
}
assert_file_empty() {
  [ ! -s "$1" ] || { echo "File '$1' is not empty but should be" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_contains_line() {
  grep -q -x "$2" "$1" || { echo "File '$1' does not contain line '$2'" && exit 1; }
}
assert_file_not_contains() {
  grep -q "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}

#################################################################

echo ">>> Prepare test data"

cat > example_R1.fastq <<'EOF'
@read1
ACGTACGTACGTAAAAA
+
IIIIIIIIIIIIIIIII
@read2
ACGTACGTACGTCCCCC
+
IIIIIIIIIIIIIIIII
EOF

cat > example_R2.fastq <<'EOF'
@read1
ACGTACGTACGTGGGGG
+
IIIIIIIIIIIIIIIII
@read2
ACGTACGTACGTTTTTT
+
IIIIIIIIIIIIIIIII
EOF

#################################################################

echo ">>> Testing for paired-end reads"
"$meta_executable" \
    --paired true \
    --input "example_R1.fastq;example_R2.fastq" \
    --adapter "ACG" \
    --trim_html_1 example_R1.trimmed.html \
    --trim_html_2 example_R2.trimmed.html \
    --trim_zip_1 example_R1.trimmed.zip \
    --trim_zip_2 example_R2.trimmed.zip \
    --fastq_1 example_R1.trimmed.fastq \
    --fastq_2 example_R2.trimmed.fastq \
    --trim_log_1 example_R1.trimming_report.txt \
    --trim_log_2 example_R2.trimming_report.txt

echo ">> Checking output"
assert_file_exists "example_R1.trimmed.html"
assert_file_exists "example_R2.trimmed.html"
assert_file_exists "example_R1.trimmed.zip"
assert_file_exists "example_R2.trimmed.zip"
assert_file_exists "example_R1.trimmed.fastq"
assert_file_exists "example_R2.trimmed.fastq"
assert_file_exists "example_R1.trimming_report.txt"
assert_file_exists "example_R2.trimming_report.txt"

echo ">> Check if output is empty"
assert_file_not_empty "example_R1.trimmed.html"
assert_file_not_empty "example_R2.trimmed.html"
assert_file_not_empty "example_R1.trimmed.zip"
assert_file_not_empty "example_R2.trimmed.zip"
assert_file_not_empty "example_R1.trimmed.fastq"
assert_file_not_empty "example_R2.trimmed.fastq"
assert_file_not_empty "example_R1.trimming_report.txt"
assert_file_not_empty "example_R2.trimming_report.txt"

echo ">> Check contents"
assert_file_contains_line "example_R1.trimmed.fastq" "TACGTACGTAAAAA"
assert_file_contains_line "example_R2.trimmed.fastq" "TACGTACGTGGGGG"
assert_file_contains "example_R1.trimming_report.txt" "sequences processed in total"
assert_file_contains "example_R2.trimming_report.txt" "Number of sequence pairs removed because at least one read was shorter than the length cutoff"

#################################################################

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --paired false \
    --input "example_R1.fastq" \
    --adapter "ACG" \
    --trim_html_1 example.trimmed.html \
    --trim_zip_1 example.trimmed.zip \
    --fastq_1 example.trimmed.fastq \
    --trim_log_1 example.trimming_report.txt \

echo ">> Checking output"
assert_file_exists "example.trimmed.html"
assert_file_exists "example.trimmed.zip"
assert_file_exists "example.trimmed.fastq"
assert_file_exists "example.trimming_report.txt"

echo ">> Check if output is empty"
assert_file_not_empty "example.trimmed.html"
assert_file_not_empty "example.trimmed.zip"
assert_file_not_empty "example.trimmed.fastq"
assert_file_not_empty "example.trimming_report.txt"

echo ">> Check contents"
assert_file_contains_line "example.trimmed.fastq" "TACGTACGTAAAAA"
assert_file_contains "example.trimming_report.txt" "Sequences removed because they became shorter than the length cutoff"

#################################################################

echo ">>> Test finished successfully"
exit 0
