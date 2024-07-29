#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/bedtools/bedtools_sort/bedtools_sort"
meta_resources_dir="src/bedtools/bedtools_sort"
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
mkdir -p test_data

# Create and populate example files
printf "chr1\t300\t400\nchr1\t150\t250\nchr1\t100\t200" > "test_data/featureA.bed"
printf "chr2\t290\t390\nchr2\t180\t280\nchr1\t500\t600" > "test_data/featureB.bed"
printf "chr3\t120\t220\nchr1\t250\t350\nchr2\t500\t580" > "test_data/featureC.bed"

# Create and populate example.gff file
printf "##gff-version 3\n" > "test_data/example.gff"
printf "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/example.gff"
printf "chr3\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/example.gff"
printf "chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=transcript1\n" >> "test_data/example.gff"
printf "chr2\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/example.gff"
printf "chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tID=cds1;Parent=transcript1\n" >> "test_data/example.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/example.gff"

# Create expected output files


# Test 1: Default sort on BED file
mkdir test1
cd test1

echo "> Run bedtools_sort on BED file"
"$meta_executable" \
  --input "../test_data/featureA.bed" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_sorted.bed"
echo "- test1 succeeded -"

cd ..

# Test 2: Default sort on GFF file
mkdir test2
cd test2

echo "> Run bedtools_sort on GFF file"
"$meta_executable" \
  --input "../test_data/example.gff" \
  --output "output.gff"

# checks
assert_file_exists "output.gff"
assert_file_not_empty "output.gff"
assert_identical_content "output.gff" "../test_data/expected_sorted.gff"
echo "- test2 succeeded -"

cd ..
