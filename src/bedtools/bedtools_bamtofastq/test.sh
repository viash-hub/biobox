#!/bin/bash

# exit on error
set -e

test_data="$meta_resources_dir/test_data"

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

# Test 1: normal conversion 
mkdir test1
cd test1

echo "> Run bedtools bamtofastq on BAM file"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --fastq "output.fastq"

# checks
assert_file_exists "output.fastq"
assert_file_not_empty "output.fastq"
assert_identical_content "output.fastq" "$test_data/expected.fastq"
echo "- test1 succeeded -"

cd ..

# Test 2: with tags
mkdir test2
cd test2

echo "> Run bedtools bamtofastq on BAM file with tags"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --fastq "output.fastq" \
  --tags

# checks
assert_file_exists "output.fastq"
assert_file_not_empty "output.fastq"
assert_identical_content "output.fastq" "$test_data/expected.fastq"
echo "- test2 succeeded -"

cd ..

# Test 3: with option fq2
mkdir test3
cd test3

echo "> Run bedtools bamtofastq on BAM file with output_fq2"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --fastq "output1.fastq" \
  --fastq2 "output2.fastq" 

# checks
assert_file_exists "output1.fastq"
assert_file_not_empty "output1.fastq"
assert_identical_content "output1.fastq" "$test_data/expected_1.fastq"
assert_file_exists "output2.fastq"
assert_file_not_empty "output2.fastq"
assert_identical_content "output2.fastq" "$test_data/expected_2.fastq"
echo "- test3 succeeded -"

cd ..

echo "All tests succeeded"
exit 0


