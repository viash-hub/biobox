#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/bedtools/bedtools_sort/bedtools_merge"
meta_resources_dir="src/bedtools/bedtools_merge"
## VIASH END

# directory of the bam file
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

# Create directories for tests
echo "Creating Test Data..."
mkdir -p test_data

# Create and populate example files
printf "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "test_data/featureA.bed"
printf "chr1\t100\t200\ta1\t1\t+\nchr1\t180\t250\ta2\t2\t+\nchr1\t250\t500\ta3\t3\t-\nchr1\t501\t1000\ta4\t4\t+\n" > "test_data/featureB.bed"

# Create and populate example.gff file
printf "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/expected_sorted.gff"
printf "chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=transcript1\n" >> "test_data/expected_sorted.gff"
printf "chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tID=cds1;Parent=transcript1\n" >> "test_data/expected_sorted.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/expected_sorted.gff"
printf "chr2\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/expected_sorted.gff"
printf "chr3\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/expected_sorted.gff"

# Create expected output files
printf "chr1\t100\t250\nchr1\t300\t400\n" > "test_data/expected.bed"
printf "chr1\t100\t250\nchr1\t250\t500\nchr1\t501\t1000\n" > "test_data/expected_strand.bed"
printf "chr1\t100\t250\nchr1\t501\t1000\n" > "test_data/expected_specific_strand.bed"
printf "chr1\t128\t228\nchr1\t428\t528\n" > "test_data/expected_bam.bed"
printf "chr1\t100\t400\n" > "test_data/expected_distance.bed"

# Test 1: Default sort on BED file
mkdir test1
cd test1

echo "> Run bedtools_merge on BED file"
"$meta_executable" \
  --input "../test_data/featureA.bed" \
  --output "output.bed"

# # checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected.bed"
echo "- test1 succeeded -"

cd ..

# Test 2: strand option
mkdir test2
cd test2

echo "> Run bedtools_merge on BED file with strand option"
"$meta_executable" \
  --input "../test_data/featureB.bed" \
  --output "output.bed" \
  --strand

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_strand.bed"
echo "- test2 succeeded -"

cd ..

# Test 3: specific strand option
mkdir test3
cd test3

echo "> Run bedtools_merge on BED file with specific strand option"
"$meta_executable" \
  --input "../test_data/featureB.bed" \
  --output "output.bed" \
  --specific_strand "+" 

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_specific_strand.bed"
echo "- test3 succeeded -"

cd ..

# Test 4: BED option
mkdir test4
cd test4

echo "> Run bedtools_merge on BAM file with BED option"
"$meta_executable" \
  --input "$test_data/feature.bam" \
  --output "output.bed" \
  --bed

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_bam.bed"
echo "- test4 succeeded -"

cd ..

# Test 5: distance option
mkdir test5
cd test5

echo "> Run bedtools_merge on BED file with distance option"
"$meta_executable" \
  --input "../test_data/featureA.bed" \
  --output "output.bed" \
  --distance -5 

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected.bed"
echo "- test5 succeeded -"

cd ..

# Test 6: columns option
# mkdir test6
# cd test6

# echo "> Run bedtools_merge on BED file with columns option"
# "$meta_executable" \
#   --input "../test_data/featureA.bed" \
#   --output "output.bed" \
#   --columns 2

# # checks
# assert_file_exists "output.bed"
# assert_file_not_empty "output.bed"
# assert_identical_content "output.bed" "../test_data/expected.bed"
# echo "- test6 succeeded -"

# cd ..


echo "---- All tests succeeded! ----"
exit 0
