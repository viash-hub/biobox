#!/bin/bash

# exit on error
set -eo pipefail

## VIASH START
# meta_executable="target/executable/fastqc"
# meta_resources_dir="src/fastqc"
## VIASH END

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
assert_identical_content() {
  diff -a "$2" "$1" \
    || (echo "Files are not identical!" && exit 1)
}
#############################################

# Create directories for tests
echo "Creating Test Data..."
TMPDIR=$(mktemp -d "$meta_temp_dir/XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Create and populate input.fasta
cat > "$TMPDIR/input_1.fq" <<EOL
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
EOL

cat > "$TMPDIR/input_2.fq" <<EOL
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
EOL

# Create and populate contaminants.txt
printf "contaminant_sequence1\tCACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGA\n" > "$TMPDIR/contaminants.txt"
printf "contaminant_sequence2\tGATCTTGG\n" >> "$TMPDIR/contaminants.txt"

# Create and populate SAM file 
printf "@HD\tVN:1.0\tSO:unsorted\n" > "$TMPDIR/example.sam"
printf "@SQ\tSN:chr1\tLN:248956422\n" >> "$TMPDIR/example.sam"
printf "@SQ\tSN:chr2\tLN:242193529\n" >> "$TMPDIR/example.sam"
printf "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.3.4.1\tCL:\"/usr/bin/bowtie2-align-s --wrapper basic-0 -x genome -U reads.fq -S output.sam\"\n" >> "$TMPDIR/example.sam"
printf "read1\t0\tchr1\t100\t255\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:-10\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:50\tYT:Z:UU\n" >> "$TMPDIR/example.sam"
printf "read2\t0\tchr2\t150\t255\t50M\t*\t0\t0\tTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:-8\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:50\tYT:Z:UU\n" >> "$TMPDIR/example.sam"
printf "read3\t16\tchr1\t200\t255\t50M\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:-12\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:50\tYT:Z:UU" >> "$TMPDIR/example.sam"

cat > "$TMPDIR/expected_summary.txt" <<EOL
PASS	Basic Statistics	input_1.fq
PASS	Per base sequence quality	input_1.fq
FAIL	Per sequence quality scores	input_1.fq
FAIL	Per base sequence content	input_1.fq
FAIL	Per sequence GC content	input_1.fq
FAIL	Per base N content	input_1.fq
PASS	Sequence Length Distribution	input_1.fq
PASS	Sequence Duplication Levels	input_1.fq
FAIL	Overrepresented sequences	input_1.fq
PASS	Adapter Content	input_1.fq
EOL

cat > "$TMPDIR/expected_summary2.txt" <<EOL
PASS	Basic Statistics	input_2.fq
PASS	Per base sequence quality	input_2.fq
FAIL	Per sequence quality scores	input_2.fq
FAIL	Per base sequence content	input_2.fq
FAIL	Per sequence GC content	input_2.fq
FAIL	Per base N content	input_2.fq
PASS	Sequence Length Distribution	input_2.fq
PASS	Sequence Duplication Levels	input_2.fq
FAIL	Overrepresented sequences	input_2.fq
PASS	Adapter Content	input_2.fq
EOL

cat > "$TMPDIR/expected_summary_sam.txt" <<EOL
PASS	Basic Statistics	example.sam
PASS	Per base sequence quality	example.sam
FAIL	Per sequence quality scores	example.sam
FAIL	Per base sequence content	example.sam
WARN	Per sequence GC content	example.sam
PASS	Per base N content	example.sam
WARN	Sequence Length Distribution	example.sam
PASS	Sequence Duplication Levels	example.sam
FAIL	Overrepresented sequences	example.sam
PASS	Adapter Content	example.sam
EOL

# Test 1: Run fastqc with default parameters
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "-> Run Test: one input"
"$meta_executable" \
    --extract \
    --input "../input_1.fq" \

assert_file_exists "input_1_fastqc.html"
assert_file_exists "input_1_fastqc.zip"
assert_file_exists "input_1_fastqc/summary.txt"
assert_file_not_empty "input_1_fastqc.html"
assert_file_not_empty "input_1_fastqc.zip"
assert_identical_content "input_1_fastqc/summary.txt" "../expected_summary.txt"
echo "- test succeeded -"

popd > /dev/null

# Test 2: Run fastqc with multiple inputs
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "-> Run Test: two inputs"
"$meta_executable" \
 --extract \
 --input "../input_1.fq" \
 --input "../input_2.fq" 

# File 1
assert_file_exists "input_1_fastqc.html"
assert_file_exists "input_1_fastqc.zip"
assert_file_exists "input_1_fastqc/summary.txt"
assert_file_not_empty "input_1_fastqc.html"
assert_file_not_empty "input_1_fastqc.zip"
assert_identical_content "input_1_fastqc/summary.txt" "../expected_summary.txt"
# File 2
assert_file_exists "input_2_fastqc.html"
assert_file_exists "input_2_fastqc.zip"
assert_file_exists "input_2_fastqc/summary.txt"
assert_file_not_empty "input_2_fastqc.html"
assert_file_not_empty "input_2_fastqc.zip"
assert_identical_content "input_2_fastqc/summary.txt" "../expected_summary2.txt"
echo "- test succeeded -"

popd > /dev/null

# Test 3: Run fastqc with contaminants
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null

echo "-> Run Test: contaminants"
"$meta_executable" \
 --extract \
 --input "../input_1.fq" \
 --contaminants "../contaminants.txt"

assert_file_exists "input_1_fastqc.html"
assert_file_exists "input_1_fastqc.zip"
assert_file_exists "input_1_fastqc/summary.txt"
assert_file_not_empty "input_1_fastqc.html"
assert_file_not_empty "input_1_fastqc.zip"
assert_identical_content "input_1_fastqc/summary.txt" "../expected_summary.txt"
assert_file_contains "contaminant" "input_1_fastqc/fastqc_data.txt"
echo "- test succeeded -"

popd > /dev/null

# Test 4: Run fastqc with sam file
mkdir "$TMPDIR/test4" && pushd "$TMPDIR/test4" > /dev/null

echo "-> Run Test: sam file"
"$meta_executable" \
 --extract \
 --input "../example.sam" \
 --format "sam"

assert_file_exists "example_fastqc.html"
assert_file_exists "example_fastqc.zip"
assert_file_exists "example_fastqc/summary.txt"
assert_file_not_empty "example_fastqc.html"
assert_file_not_empty "example_fastqc.zip"
assert_identical_content "example_fastqc/summary.txt" "../expected_summary_sam.txt"
echo "- test succeeded -"

popd > /dev/null

# Test 5: Run fastqc with multiple options
mkdir "$TMPDIR/test5" && pushd "$TMPDIR/test5" > /dev/null

echo "-> Run Test: multiple options"
"$meta_executable" \
 --extract \
 --input "../input_1.fq" \
 --contaminants "../contaminants.txt" \
 --format "fastq" \
 --casava \
 --nofilter \
 --nogroup \
 --min_length 10 \
 --kmers 5

assert_file_exists "input_1_fastqc.html"
assert_file_exists "input_1_fastqc.zip"
assert_file_exists "input_1_fastqc/summary.txt"
assert_file_not_empty "input_1_fastqc.html"
assert_file_not_empty "input_1_fastqc.zip"
assert_identical_content "input_1_fastqc/summary.txt" "../expected_summary.txt"
assert_file_contains "contaminant" "input_1_fastqc/fastqc_data.txt"
echo "- test succeeded -"

popd > /dev/null

# Test 6: Run fastqc with output options
mkdir "$TMPDIR/test6" && pushd "$TMPDIR/test6" > /dev/null

echo "-> Run Test: output options"
"$meta_executable" \
 --input "../input_1.fq" \
 --html "../output_*_.html" \
 --zip "../output_*_.zip"

# Check if the html file was generated
[ ! -f "test_data/output_input_1_.html" ] \
    && echo "Output HTML file not found." && exit 1

# Check if the zip file was generated
[ ! -f "test_data/output_input_1_.zip" ] \
    && echo "Output ZIP file not found." && exit 1

# Check if the files are empty
[ ! -s "test_data/output_input_1_.html" ] \
    && echo "Output HTML file is empty." && exit 1

[ ! -s "test_data/output_input_1_.zip" ] \
    && echo "Output ZIP file is empty." && exit 1

assert_file_exists "output_1_fastqc.html"
assert_file_exists "output_1_fastqc.zip"
assert_file_not_empty "output_1_fastqc.html"
assert_file_not_empty "output_1_fastqc.zip"
echo "- test succeeded -"

popd > /dev/null

echo "All tests succeeded!"
exit 0