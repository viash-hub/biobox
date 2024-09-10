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
assert_file_not_contains() {
  grep -q "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}

#################################################################

echo ">>> Prepare test data"

cat > example_R1.fastq <<'EOF'
@SRR6357071.22842410 22842410/1 kraken:taxid|4932
CAAGTTTTCATCTTCAACAGCTGATTGACTTCTTTGTGGTATGCCTCGATATATTTTTCTTTTTCTTTAATATCTTTATTATAGGTGATTGCCTCATCGTA
+
BBBBBFFFFFFFFFFFFFFF/BFFFFFFFFFFFFFFFFBFFBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFBF<
@SRR6357071.52260105 52260105/1 kraken:taxid|4932
TAGACTTACCAGTACCCTTTTCGACGGCGGAAACATTCAAAATACCGTTAGAGTCGACATCGAAAGTGACTTCAATTTGTGGGACACCTCTTGGAGCTGGT
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/FFFFFFFFFFFFFFFF
EOF

cat > example_R2.fastq <<'EOF'
@SRR6357071.22842410 22842410/2 kraken:taxid|4932
CCGAGATCGAAGAAACGAATTCACCTGATTGCAGCTGTAAAAGCAGTAAAATCAATCAAACCAATACGGACAACCTTACGATACGATGAGGCAATCACCTA
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SRR6357071.52260105 52260105/2 kraken:taxid|4932
GTTGATTCCAAGAAACTCTACCATTCCAACTAAGAAATCCGAAGTTTTCTCTACTTATGCTGACAACCAACCAGGTGTCTTGATTCAAGTCTTTGAAGGTG
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
EOF

#################################################################

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --input "example_R1.fastq" \
    --trimmed_fastqc_html_1 output_se_test/example.trimmed.html \
    --trimmed_fastqc_zip_1 output_se_test/example.trimmed.zip \
    --trimmed_r1 output_se_test/example.trimmed.fastq \
    --trimming_report_r1 output_se_test/example.trimming_report.txt \
    --fastqc \
    --output_dir output_se_test

echo ">> Checking output"
assert_file_exists "output_se_test/example.trimmed.html"
assert_file_exists "output_se_test/example.trimmed.zip"
assert_file_exists "output_se_test/example.trimmed.fastq"
assert_file_exists "output_se_test/example.trimming_report.txt"

echo ">> Check if output is empty"
assert_file_not_empty "output_se_test/example.trimmed.html"
assert_file_not_empty "output_se_test/example.trimmed.zip"
assert_file_not_empty "output_se_test/example.trimmed.fastq"
assert_file_not_empty "output_se_test/example.trimming_report.txt"

echo ">> Check contents"
assert_file_contains "output_se_test/example.trimmed.fastq" "@SRR6357071.22842410 22842410/1"
assert_file_contains "output_se_test/example.trimming_report.txt" "Sequences removed because they became shorter than the length cutoff"

#################################################################

echo ">>> Testing for paired-end reads"
"$meta_executable" \
    --paired \
    --input "example_R1.fastq;example_R2.fastq" \
    --trimmed_fastqc_html_1 output_pe_test/example_R1.trimmed.html \
    --trimmed_fastqc_html_2 output_pe_test/example_R2.trimmed.html \
    --trimmed_fastqc_zip_1 output_pe_test/example_R1.trimmed.zip \
    --trimmed_fastqc_zip_2 output_pe_test/example_R2.trimmed.zip \
    --trimmed_r1 output_pe_test/example_R1.trimmed.fastq \
    --trimmed_r2 output_pe_test/example_R2.trimmed.fastq \
    --trimming_report_r1 output_pe_test/example_R1.trimming_report.txt \
    --trimming_report_r2 output_pe_test/example_R2.trimming_report.txt \
    --fastqc \
    --output_dir output_pe_test

echo ">> Checking output"
assert_file_exists "output_pe_test/example_R1.trimmed.html"
assert_file_exists "output_pe_test/example_R2.trimmed.html"
assert_file_exists "output_pe_test/example_R1.trimmed.zip"
assert_file_exists "output_pe_test/example_R2.trimmed.zip"
assert_file_exists "output_pe_test/example_R1.trimmed.fastq"
assert_file_exists "output_pe_test/example_R2.trimmed.fastq"
assert_file_exists "output_pe_test/example_R1.trimming_report.txt"
assert_file_exists "output_pe_test/example_R2.trimming_report.txt"

echo ">> Check if output is empty"
assert_file_not_empty "output_pe_test/example_R1.trimmed.html"
assert_file_not_empty "output_pe_test/example_R2.trimmed.html"
assert_file_not_empty "output_pe_test/example_R1.trimmed.zip"
assert_file_not_empty "output_pe_test/example_R2.trimmed.zip"
assert_file_not_empty "output_pe_test/example_R1.trimmed.fastq"
assert_file_not_empty "output_pe_test/example_R2.trimmed.fastq"
assert_file_not_empty "output_pe_test/example_R1.trimming_report.txt"
assert_file_not_empty "output_pe_test/example_R2.trimming_report.txt"

echo ">> Check contents"
assert_file_contains "output_pe_test/example_R1.trimmed.fastq" "@SRR6357071.22842410 22842410/1"
assert_file_contains "output_pe_test/example_R2.trimmed.fastq" "@SRR6357071.22842410 22842410/2"
assert_file_contains "output_pe_test/example_R1.trimming_report.txt" "sequences processed in total"
assert_file_contains "output_pe_test/example_R2.trimming_report.txt" "Number of sequence pairs removed because at least one read was shorter than the length cutoff"

#################################################################

echo ">>> Test finished successfully"
exit 0
