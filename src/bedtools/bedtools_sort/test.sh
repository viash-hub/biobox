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
printf "#Header\nchr1\t300\t400\nchr1\t150\t250\nchr1\t100\t200" > "test_data/featureA.bed"
printf "chr2\t290\t400\nchr2\t180\t220\nchr1\t500\t600" > "test_data/featureB.bed"
printf "chr1\t100\t200\tfeature1\t960\nchr1\t150\t250\tfeature2\t850\nchr1\t300\t400\tfeature3\t740\nchr2\t290\t390\tfeature4\t630\nchr2\t180\t280\tfeature5\t920\nchr3\t120\t220\tfeature6\t410\n" > "test_data/featureC.bed"
printf "chr1\nchr3\nchr2\n" > "test_data/genome.txt"
printf "chr1\t248956422\nchr3\t242193529\nchr2\t198295559\n" > "test_data/genome.fai"

# Create and populate example.gff file
printf "##gff-version 3\n" > "test_data/example.gff"
printf "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/example.gff"
printf "chr3\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/example.gff"
printf "chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=transcript1\n" >> "test_data/example.gff"
printf "chr2\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/example.gff"
printf "chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tID=cds1;Parent=transcript1\n" >> "test_data/example.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/example.gff"

# Create expected output files
printf "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "test_data/expected_sorted_A.bed"
printf "chr2\t180\t220\nchr1\t500\t600\nchr2\t290\t400\n" > "test_data/expected_sizeA.bed"
printf "chr2\t290\t400\nchr1\t500\t600\nchr2\t180\t220\n" > "test_data/expected_sizeD.bed"
printf "chr1\t500\t600\nchr2\t180\t220\nchr2\t290\t400\n" > "test_data/expected_chrThenSizeA.bed"
printf "chr1\t500\t600\nchr2\t290\t400\nchr2\t180\t220\n" > "test_data/expected_chrThenSizeD.bed"
printf "chr1\t300\t400\tfeature3\t740\nchr1\t150\t250\tfeature2\t850\nchr1\t100\t200\tfeature1\t960\nchr2\t290\t390\tfeature4\t630\nchr2\t180\t280\tfeature5\t920\nchr3\t120\t220\tfeature6\t410\n" > "test_data/expected_chrThenScoreA.bed"
printf "chr1\t100\t200\tfeature1\t960\nchr1\t150\t250\tfeature2\t850\nchr1\t300\t400\tfeature3\t740\nchr2\t180\t280\tfeature5\t920\nchr2\t290\t390\tfeature4\t630\nchr3\t120\t220\tfeature6\t410\n" > "test_data/expected_chrThenScoreD.bed"
printf "chr1\t100\t200\tfeature1\t960\nchr1\t150\t250\tfeature2\t850\nchr1\t300\t400\tfeature3\t740\nchr3\t120\t220\tfeature6\t410\nchr2\t180\t280\tfeature5\t920\nchr2\t290\t390\tfeature4\t630\n" > "test_data/expected_genome.bed"
printf "#Header\nchr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "test_data/expected_header.bed"

# expected_sorted.gff
printf "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/expected_sorted.gff"
printf "chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=transcript1\n" >> "test_data/expected_sorted.gff"
printf "chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tID=cds1;Parent=transcript1\n" >> "test_data/expected_sorted.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/expected_sorted.gff"
printf "chr2\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/expected_sorted.gff"
printf "chr3\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/expected_sorted.gff"

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
assert_identical_content "output.bed" "../test_data/expected_sorted_A.bed"
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

# Test 3: Sort on sizeA
mkdir test3
cd test3

echo "> Run bedtools_sort on BED file with sizeA"
"$meta_executable" \
  --input "../test_data/featureB.bed" \
  --output "output.bed" \
  --sizeA

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_sizeA.bed"
echo "- test3 succeeded -"

cd ..

# Test 4: Sort on sizeD
mkdir test4
cd test4

echo "> Run bedtools_sort on BED file with sizeD"
"$meta_executable" \
  --input "../test_data/featureB.bed" \
  --output "output.bed" \
  --sizeD

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_sizeD.bed"
echo "- test4 succeeded -"

cd ..

# Test 5: Sort on chrThenSizeA
mkdir test5
cd test5

echo "> Run bedtools_sort on BED file with chrThenSizeA"
"$meta_executable" \
    --input "../test_data/featureB.bed" \
    --output "output.bed" \
    --chrThenSizeA

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_chrThenSizeA.bed"
echo "- test5 succeeded -"

cd ..

# Test 6: Sort on chrThenSizeD
mkdir test6
cd test6

echo "> Run bedtools_sort on BED file with chrThenSizeD"
"$meta_executable" \
    --input "../test_data/featureB.bed" \
    --output "output.bed" \
    --chrThenSizeD

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_chrThenSizeD.bed"
echo "- test6 succeeded -"

cd ..

# Test 7: Sort on chrThenScoreA
mkdir test7
cd test7

echo "> Run bedtools_sort on BED file with chrThenScoreA"
"$meta_executable" \
    --input "../test_data/featureC.bed" \
    --output "output.bed" \
    --chrThenScoreA

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_chrThenScoreA.bed"
echo "- test7 succeeded -"

cd ..

# Test 8: Sort on chrThenScoreD
mkdir test8
cd test8

echo "> Run bedtools_sort on BED file with chrThenScoreD"
"$meta_executable" \
    --input "../test_data/featureC.bed" \
    --output "output.bed" \
    --chrThenScoreD

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_chrThenScoreD.bed"
echo "- test8 succeeded -"

cd ..

# Test 9: Sort according to genome file
mkdir test9
cd test9

echo "> Run bedtools_sort on BED file according to genome file"
"$meta_executable" \
    --input "../test_data/featureC.bed" \
    --output "output.bed" \
    --genome "../test_data/genome.txt"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_genome.bed"
echo "- test9 succeeded -"

cd ..

# Test 10: Sort according to faidx file
mkdir test10
cd test10

echo "> Run bedtools_sort on BED file according to faidx file"
"$meta_executable" \
    --input "../test_data/featureC.bed" \
    --output "output.bed" \
    --faidx "../test_data/genome.fai"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_genome.bed"
echo "- test10 succeeded -"

cd ..

# Test 11: Sort with header
mkdir test11
cd test11

echo "> Run bedtools_sort on BED file with header"
"$meta_executable" \
    --input "../test_data/featureA.bed" \
    --output "output.bed" \
    --header

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_header.bed"
echo "- test11 succeeded -"

cd ..

echo "---- All tests succeeded! ----"
exit 0
