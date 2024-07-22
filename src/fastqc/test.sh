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

# Test 1: Run fastqc with default parameters
mkdir test1
cd test1
echo "-> Run Test without options"
"$meta_executable" \
 --input "../test_data/input_1.fq" 

# Check if the html file was generated
[ ! -f "../test_data/input_1_fastqc.html" ] \
    && echo "Output HTML file not found." && exit 1

# Check if the zip file was generated
[ ! -f "../test_data/input_1_fastqc.zip" ] \
    && echo "Output ZIP file not found." && exit 1

# Check if the files are empty
[ ! -s "../test_data/input_1_fastqc.html" ] \
    && echo "Output HTML file is empty." && exit 1

[ ! -s "../test_data/input_1_fastqc.zip" ] \
    && echo "Output ZIP file is empty." && exit 1

# Unzip the zip file
unzip "../test_data/input_1_fastqc.zip" -d "../test_data/test1/"

# Check if the summary.txt was extracted
[ ! -f "../test_data/test1/input_1_fastqc/summary.txt" ] && echo "Extracted files not found." && exit 1

# Check if the summary.txt is correct
diff -a "../test_data/expected_summary.txt" "../test_data/test1/input_1_fastqc/summary.txt" \
    || (echo "Output summary file does not match expected output" && exit 1)

echo "- test succeeded -"
cd ..

# Test 2: Run fastqc with multiple inputs
mkdir test2
cd test2
echo "-> Run Test with multiple inputs"
"$meta_executable" \
 --input "../test_data/input_1.fq,../test_data/input_2.fq"

# Check if the html files was generated
[ ! -f "../test_data/input_1_fastqc.html" ] && [ ! -f "../test_data/input_2_fastqc.html" ] \
    && echo "Output HTML files not found." && exit 1

# Check if the zip files was generated
[ ! -f "../test_data/input_1_fastqc.zip" ] && [ ! -f "../test_data/input_2_fastqc.zip" ] \
    && echo "Output ZIP files not found." && exit 1

# Check if the files are empty
[ ! -s "../test_data/input_1_fastqc.html" ] && [ ! -s "../test_data/input_2_fastqc.html" ] \
    && echo "Output HTML files are empty." && exit 1

[ ! -s "../test_data/input_1_fastqc.zip" ] && [ ! -s "../test_data/input_2_fastqc.zip" ] \
    && echo "Output ZIP files are empty." && exit 1

# Unzip the zip files
unzip "../test_data/input_1_fastqc.zip" -d "../test_data/test2/"
unzip "../test_data/input_2_fastqc.zip" -d "../test_data/test2/"

# Check if the summary.txt was extracted
[ ! -f "../test_data/test2/input_1_fastqc/summary.txt" ] && echo "Extracted files not found." && exit 1
[ ! -f "../test_data/test2/input_2_fastqc/summary.txt" ] && echo "Extracted files not found." && exit 1

# Check if the summary.txt is correct
diff -a "../test_data/expected_summary.txt" "../test_data/test2/input_1_fastqc/summary.txt" \
    || (echo "Output summary file does not match expected output" && exit 1)
diff -a "../test_data/expected_summary2.txt" "../test_data/test2/input_2_fastqc/summary.txt" \
    || (echo "Output summary file does not match expected output" && exit 1)

echo "- test succeeded -"
cd ..

echo "All tests succeeded!"
exit 0