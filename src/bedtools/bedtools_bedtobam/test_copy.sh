#!/bin/bash

# exit on error
set -eo pipefail

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
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\t128\t228\t255,0,0\t1\t100\t0\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\t428\t528\t255,0,0\t1\t100\t0\n" > "test_data/example.bed12"
# Create and populate example.gff file
printf "##gff-version 3\n" > "test_data/example.gff"
printf "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/example.gff"
printf "chr3\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/example.gff"
printf "chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=transcript1\n" >> "test_data/example.gff"
printf "chr2\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/example.gff"
printf "chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tID=cds1;Parent=transcript1\n" >> "test_data/example.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/example.gff"

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
echo "- test1 succeeded -"

cd ..

# I could not assert_identical_content (bam file is compressed and -ubam option does not seem to work on this version).
# could use samtools view to check the content of the bam file, would need to add samtools as a dependency.
# TODO: check if older versions works with -ubam option!!!
# Version 2.27.1 (bedtools) -ubam outputs better than version 2.30.0 (bedtools).
# However, its output is still somewhat compressed.

# Test 2: BED12 file
mkdir test2
cd test2

echo "> Run bedtools_bedtobam on BED12 file"
"$meta_executable" \
  --input "../test_data/example.bed12" \
  --genome "../test_data/genome.txt" \
  --output "output.bam" \
  --bed12 \
#  --uncompress_bam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
# assert_identical_content "output.bam" "../test_data/expected12.sam"
echo "- test2 succeeded -"

cd ..

# Test 3: Uncompressed BAM file
mkdir test3
cd test3

echo "> Run bedtools_bedtobam on BED file with uncompressed BAM output"
"$meta_executable" \
  --input "../test_data/example.bed" \
  --genome "../test_data/genome.txt" \
  --output "output.bam" \
  --uncompress_bam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
# assert_identical_content "output.bam" "../test_data/expected.sam"

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
  --map_quality 10 \
#  --uncompress_bam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
#assert_identical_content "output.bam" "../test_data/expected_mapquality.sam"
echo "- test4 succeeded -"

cd ..

# Test 5: gff to bam conversion
mkdir test5
cd test5

echo "> Run bedtools_bedtobam on GFF file"
"$meta_executable" \
  --input "../test_data/example.gff" \
  --genome "../test_data/genome.txt" \
  --output "output.bam"
#  --uncompress_bam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
# assert_identical_content "output.bam" "../test_data/expected_gff.sam"
echo "- test5 succeeded -"

cd ..

echo "---- All tests succeeded! ----"
exit 0
