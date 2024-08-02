#!/bin/bash

# exit on error
set -e

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
mkdir -p test_data

# Create and populate input files
printf "chr1\t248956422\nchr3\t242193529\nchr2\t198295559\n" > "test_data/genome.txt"
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t37\t+\nchr2:172936693-172938111\t428\t528\tmy_read/2\t37\t-\n" > "test_data/example.bed"
printf "" > "test_data/example.bed12"

# Test 1: Default conversion BED to BAM
mkdir test1
cd test1

echo "> Run bedtools_bedtobam on BED file"
"$meta_executable" \
  --input "../test_data/example.bed" \
  --genome "../test_data/genome.txt" \
  --output "output.bam"

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
# assert_identical_content # (bam file is compressed and ubam option does not seem to work.)
# could use samtools view to check the content of the bam file, would need to add samtools as a dependency
echo "- test1 succeeded -"

cd ..

# Test 2: BED12 file
mkdir test2
cd test2

echo "> Run bedtools_bedtobam on BED12 file"
"$meta_executable" \
  --input "../test_data/example.bed12" \
  --genome "../test_data/genome.txt" \
  --output "output.bam" \
  --bed12

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"

echo "- test2 succeeded -"

cd ..

# Test 3: Uncompressed BAM file
mkdir test3
cd test3

# TODO: check if older versions works with -ubam option!!!

echo "> Run bedtools_bedtobam on BED file with uncompressed BAM output"
"$meta_executable" \
  --input "../test_data/example.bed" \
  --genome "../test_data/genome.txt" \
  --output "output.bam" \
  --uncompress_bam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"

echo "- test3 succeeded -"

cd ..

# Test 4: Map quality
mkdir test4
cd test4

echo "> Run bedtools_bedtobam on BED file with map quality"
"$meta_executable" \
  --input "../test_data/example.bed" \
  --genome "../test_data/genome.txt" \
  --output "output.bam" \
  --map_quality 10

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"

echo "- test4 succeeded -"

cd ..

echo "---- All tests succeeded! ----"
exit 0
