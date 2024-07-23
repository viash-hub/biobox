#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/fastqc"
meta_resources_dir="src/fastqc"
## VIASH END

# Create directories for tests
echo "Creating Test Data..."
mkdir test_data

# Create and populate input.fasta
cat > "test_data/input_1.fq" <<EOL
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
EOL

cat > "test_data/input_2.fq" <<EOL
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
EOL

# Create and populate contaminants.txt
printf "contaminant_sequence1\tCACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGA\n" > "test_data/contaminants.txt"
printf "contaminant_sequence2\tGATCTTGG\n" >> "test_data/contaminants.txt"

# Create and populate SAM file 
printf "@HD\tVN:1.0\tSO:unsorted\n" > "test_data/example.sam"
printf "@SQ\tSN:chr1\tLN:248956422\n" >> "test_data/example.sam"
printf "@SQ\tSN:chr2\tLN:242193529\n" >> "test_data/example.sam"
printf "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.3.4.1\tCL:\"/usr/bin/bowtie2-align-s --wrapper basic-0 -x genome -U reads.fq -S output.sam\"\n" >> "test_data/example.sam"
printf "read1\t0\tchr1\t100\t255\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:-10\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:50\tYT:Z:UU\n" >> "test_data/example.sam"
printf "read2\t0\tchr2\t150\t255\t50M\t*\t0\t0\tTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:-8\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:50\tYT:Z:UU\n" >> "test_data/example.sam"
printf "read3\t16\tchr1\t200\t255\t50M\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:-12\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:50\tYT:Z:UU" >> "test_data/example.sam"

cat > "test_data/expected_summary.txt" <<EOL
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

cat > "test_data/expected_summary2.txt" <<EOL
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

cat > "test_data/expected_summary_sam.txt" <<EOL
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
echo "-> Run Test: one input"
"$meta_executable" \
 --extract \
 --input "test_data/input_1.fq" 

# Check if the html file was generated
[ ! -f "test_data/input_1_fastqc.html" ] \
    && echo "Output HTML file not found." && exit 1

# Check if the zip file was generated
[ ! -f "test_data/input_1_fastqc.zip" ] \
    && echo "Output ZIP file not found." && exit 1

# Check if the files are empty
[ ! -s "test_data/input_1_fastqc.html" ] \
    && echo "Output HTML file is empty." && exit 1

[ ! -s "test_data/input_1_fastqc.zip" ] \
    && echo "Output ZIP file is empty." && exit 1

# Check if the summary.txt was extracted
[ ! -f "test_data/input_1_fastqc/summary.txt" ] && echo "Extracted files not found." && exit 1

# Check if the summary.txt is correct
diff -a "test_data/expected_summary.txt" "test_data/input_1_fastqc/summary.txt" \
    || (echo "Output summary file does not match expected output" && exit 1)

rm -r "test_data/input_1_fastqc"
rm "test_data/input_1_fastqc.html"
rm "test_data/input_1_fastqc.zip"

echo "- test succeeded -"

# Test 2: Run fastqc with multiple inputs
echo "-> Run Test: two inputs"
"$meta_executable" \
 --extract \
 --input "test_data/input_1.fq,test_data/input_2.fq"

# Check if the html files was generated
[ ! -f "test_data/input_1_fastqc.html" ] && [ ! -f "test_data/input_2_fastqc.html" ] \
    && echo "Output HTML files not found." && exit 1

# Check if the zip files was generated
[ ! -f "test_data/input_1_fastqc.zip" ] && [ ! -f "test_data/input_2_fastqc.zip" ] \
    && echo "Output ZIP files not found." && exit 1

# Check if the files are empty
[ ! -s "test_data/input_1_fastqc.html" ] && [ ! -s "test_data/input_2_fastqc.html" ] \
    && echo "Output HTML files are empty." && exit 1

[ ! -s "test_data/input_1_fastqc.zip" ] && [ ! -s "test_data/input_2_fastqc.zip" ] \
    && echo "Output ZIP files are empty." && exit 1

# Check if the summary.txt was extracted
[ ! -f "test_data/input_1_fastqc/summary.txt" ] && echo "Extracted files not found." && exit 1
[ ! -f "test_data/input_2_fastqc/summary.txt" ] && echo "Extracted files not found." && exit 1

# Check if the summary.txt is correct
diff -a "test_data/expected_summary.txt" "test_data/input_1_fastqc/summary.txt" \
    || (echo "Output summary file does not match expected output" && exit 1)
diff -a "test_data/expected_summary2.txt" "test_data/input_2_fastqc/summary.txt" \
    || (echo "Output summary file does not match expected output" && exit 1)

rm -r "test_data/input_1_fastqc"
rm -r "test_data/input_2_fastqc"
rm "test_data/input_1_fastqc.html"
rm "test_data/input_2_fastqc.html"
rm "test_data/input_1_fastqc.zip"
rm "test_data/input_2_fastqc.zip"

echo "- test succeeded -"

# Test 3: Run fastqc with contaminants
echo "-> Run Test: contaminants"
"$meta_executable" \
 --extract \
 --input "test_data/input_1.fq" \
 --contaminants "test_data/contaminants.txt"

# Check if the html file was generated
[ ! -f "test_data/input_1_fastqc.html" ] \
    && echo "Output HTML file not found." && exit 1

# Check if the zip file was generated
[ ! -f "test_data/input_1_fastqc.zip" ] \
    && echo "Output ZIP file not found." && exit 1

# Check if the files are empty
[ ! -s "test_data/input_1_fastqc.html" ] \
    && echo "Output HTML file is empty." && exit 1

[ ! -s "test_data/input_1_fastqc.zip" ] \
    && echo "Output ZIP file is empty." && exit 1

# Check if the summary.txt was extracted
[ ! -f "test_data/input_1_fastqc/summary.txt" ] && echo "Extracted files not found." && exit 1

# Checking for contaminants in fastqc_data.txt
echo "Checking for contaminants in fastqc_data.txt"
result=$(cat test_data/input_1_fastqc/fastqc_data.txt | grep "contaminant" )
expecte_result=$(printf "CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNN\t1\t100.0\tcontaminant_sequence1 (100%% over 44bp)\n")

[ -z "$result" ] && echo "Contaminants not found in fastqc_data.txt" && exit 1

[ "$result" != "$expecte_result" ] \
 && echo "Contaminants do not match expected output" \
 && echo "Result: $result" \
 && echo "Expected: $expecte_result" \
 && exit 1

rm -r "test_data/input_1_fastqc"
rm "test_data/input_1_fastqc.html"
rm "test_data/input_1_fastqc.zip"

echo "- test succeeded -"

# Test 4: Run fastqc with sam file
echo "-> Run Test: sam file"
"$meta_executable" \
 --extract \
 --input "test_data/example.sam" \
 --format "sam"

# Check if the html file was generated
[ ! -f "test_data/example_fastqc.html" ] \
    && echo "Output HTML file not found." && exit 1

# Check if the zip file was generated
[ ! -f "test_data/example_fastqc.zip" ] \
    && echo "Output ZIP file not found." && exit 1

# Check if the files are empty
[ ! -s "test_data/example_fastqc.html" ] \
    && echo "Output HTML file is empty." && exit 1

[ ! -s "test_data/example_fastqc.zip" ] \
    && echo "Output ZIP file is empty." && exit 1

# Check if the summary.txt was extracted
[ ! -f "test_data/example_fastqc/summary.txt" ] && echo "Extracted files not found." && exit 1

# Check if the summary.txt is correct
diff -a "test_data/expected_summary_sam.txt" "test_data/example_fastqc/summary.txt" \
    || (echo "Output summary file does not match expected output" && exit 1)

rm -r "test_data/example_fastqc"
rm "test_data/example_fastqc.html"
rm "test_data/example_fastqc.zip"

echo "- test succeeded -"

# Test 5: Run fastqc with multiple options
echo "-> Run Test: multiple options"
"$meta_executable" \
 --extract \
 --input "test_data/input_1.fq" \
 --contaminants "test_data/contaminants.txt" \
 --format "fastq" \
 --casava \
 --nofilter \
 --nogroup \
 --min_length 10 \
 --threads 4 \
 --kmers 5


# Check if the html file was generated
[ ! -f "test_data/input_1_fastqc.html" ] \
    && echo "Output HTML file not found." && exit 1

# Check if the zip file was generated
[ ! -f "test_data/input_1_fastqc.zip" ] \
    && echo "Output ZIP file not found." && exit 1

# Check if the files are empty
[ ! -s "test_data/input_1_fastqc.html" ] \
    && echo "Output HTML file is empty." && exit 1

[ ! -s "test_data/input_1_fastqc.zip" ] \
    && echo "Output ZIP file is empty." && exit 1

# Check if the summary.txt was extracted
[ ! -f "test_data/input_1_fastqc/summary.txt" ] && echo "Extracted files not found." && exit 1

# Check if the summary.txt is correct
diff -a "test_data/expected_summary.txt" "test_data/input_1_fastqc/summary.txt" \
    || (echo "Output summary file does not match expected output" && exit 1)

rm -r "test_data/input_1_fastqc"
rm "test_data/input_1_fastqc.html"
rm "test_data/input_1_fastqc.zip"

echo "- test succeeded -"


# Add more tests here

echo "All tests succeeded!"
exit 0