#!/bin/bash

set -e
set -eo pipefail

#############################################
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
#############################################

mkdir test_multiple_output
cd test_multiple_output

echo "#############################################"
echo "> Run cutadapt with multiple outputs"

cat > example.fa <<'EOF'
>read1
MYSEQUENCEADAPTER
>read2
MYSEQUENCEADAP
>read3
MYSEQUENCEADAPTERSOMETHINGELSE
>read4
MYSEQUENCEADABTER
>read5
MYSEQUENCEADAPTR
>read6
MYSEQUENCEADAPPTER
>read7
ADAPTERMYSEQUENCE
>read8
PTERMYSEQUENCE
>read9
SOMETHINGADAPTERMYSEQUENCE
EOF

"$meta_executable" \
  --report minimal \
  --output "out_test/*.fasta" \
  --adapter ADAPTER \
  --input example.fa \
  --fasta \
  --demultiplex_mode single \
  --no_match_adapter_wildcards \
  --json

echo ">> Checking output"
assert_file_exists "report.json"
assert_file_exists "out_test/1_001.fasta"
assert_file_exists "out_test/unknown_001.fasta"

cd ..
echo

#############################################
mkdir test_simple_single_end
cd test_simple_single_end

echo "#############################################"
echo "> Run cutadapt on single-end data"

cat > example.fa <<'EOF'
>read1
MYSEQUENCEADAPTER
>read2
MYSEQUENCEADAP
>read3
MYSEQUENCEADAPTERSOMETHINGELSE
>read4
MYSEQUENCEADABTER
>read5
MYSEQUENCEADAPTR
>read6
MYSEQUENCEADAPPTER
>read7
ADAPTERMYSEQUENCE
>read8
PTERMYSEQUENCE
>read9
SOMETHINGADAPTERMYSEQUENCE
EOF

"$meta_executable" \
  --report minimal \
  --output "out_test1/*.fasta" \
  --adapter ADAPTER \
  --input example.fa \
  --demultiplex_mode single \
  --fasta \
  --no_match_adapter_wildcards \
  --json

echo ">> Checking output"
assert_file_exists "report.json"
assert_file_exists "out_test1/1_001.fasta"
assert_file_exists "out_test1/unknown_001.fasta"

echo ">> Check if output is empty"
assert_file_not_empty "report.json"
assert_file_not_empty "out_test1/1_001.fasta"
assert_file_not_empty "out_test1/unknown_001.fasta"

echo ">> Check contents"
for i in 1 2 3 7 9; do
  assert_file_contains "out_test1/1_001.fasta" ">read$i"
done
for i in 4 5 6 8; do
  assert_file_contains "out_test1/unknown_001.fasta" ">read$i"
done

cd ..
echo

#############################################
mkdir test_multiple_single_end
cd test_multiple_single_end

echo "#############################################"
echo "> Run with a combination of inputs"

cat > example.fa <<'EOF'
>read1
ACGTACGTACGTAAAAA
>read2
ACGTACGTACGTCCCCC
>read3
ACGTACGTACGTGGGGG
>read4
ACGTACGTACGTTTTTT
EOF

cat > adapters1.fasta <<'EOF'
>adapter1
CCCCC
EOF

cat > adapters2.fasta <<'EOF'
>adapter2
GGGGG
EOF

"$meta_executable" \
  --report minimal \
  --output "out_test2/*.fasta" \
  --adapter AAAAA \
  --adapter_fasta adapters1.fasta \
  --adapter_fasta adapters2.fasta \
  --demultiplex_mode single \
  --input example.fa \
  --fasta \
  --json

echo ">> Checking output"
assert_file_exists "report.json"
assert_file_exists "out_test2/1_001.fasta"
assert_file_exists "out_test2/adapter1_001.fasta"
assert_file_exists "out_test2/adapter2_001.fasta"
assert_file_exists "out_test2/unknown_001.fasta"

echo ">> Check if output is empty"
assert_file_not_empty "report.json"
assert_file_not_empty "out_test2/1_001.fasta"
assert_file_not_empty "out_test2/adapter1_001.fasta"
assert_file_not_empty "out_test2/adapter2_001.fasta"
assert_file_not_empty "out_test2/unknown_001.fasta"

echo ">> Check contents"
assert_file_contains "out_test2/1_001.fasta" ">read1"
assert_file_contains "out_test2/adapter1_001.fasta" ">read2"
assert_file_contains "out_test2/adapter2_001.fasta" ">read3"
assert_file_contains "out_test2/unknown_001.fasta" ">read4"

cd ..
echo

#############################################
mkdir test_simple_paired_end
cd test_simple_paired_end

echo "#############################################"
echo "> Run cutadapt on paired-end data"

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

"$meta_executable" \
  --report minimal \
  --output "out_test3/*.fastq" \
  --adapter AAAAA \
  --adapter_r2 GGGGG \
  --input example_R1.fastq \
  --input_r2 example_R2.fastq \
  --quality_cutoff 20 \
  --demultiplex_mode unique_dual \
  --json \
  ---cpus 1

echo ">> Checking output"
assert_file_exists "report.json"
assert_file_exists "out_test3/1_R1_001.fastq"
assert_file_exists "out_test3/1_R2_001.fastq"
assert_file_exists "out_test3/unknown_R1_001.fastq"
assert_file_exists "out_test3/unknown_R2_001.fastq"

echo ">> Check if output is empty"
assert_file_not_empty "report.json"
assert_file_not_empty "out_test3/1_R1_001.fastq"
assert_file_not_empty "out_test3/1_R2_001.fastq"
assert_file_not_empty "out_test3/unknown_R1_001.fastq"

echo ">> Check contents"
assert_file_contains "out_test3/1_R1_001.fastq" "@read1"
assert_file_contains "out_test3/1_R2_001.fastq" "@read1"
assert_file_contains "out_test3/unknown_R1_001.fastq" "@read2"
assert_file_contains "out_test3/unknown_R2_001.fastq" "@read2"

cd ..
echo

#############################################

echo "#############################################"
echo "> Test successful"

