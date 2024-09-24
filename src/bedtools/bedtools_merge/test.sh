#!/bin/bash

# exit on error
set -eo pipefail

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
TMPDIR=$(mktemp -d "$meta_temp_dir/XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Create and populate example files
printf "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "$TMPDIR/featureA.bed"
printf "chr1\t100\t200\ta1\t1\t+\nchr1\t180\t250\ta2\t2\t+\nchr1\t250\t500\ta3\t3\t-\nchr1\t501\t1000\ta4\t4\t+\n" > "$TMPDIR/featureB.bed"
printf "chr1\t100\t200\ta1\t1.9\t+\nchr1\t180\t250\ta2\t2.5\t+\nchr1\t250\t500\ta3\t3.3\t-\nchr1\t501\t1000\ta4\t4\t+\n" > "$TMPDIR/feature_precision.bed"

# Create and populate feature.gff file
printf "##gff-version 3\n" > "$TMPDIR/feature.gff"
printf "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "$TMPDIR/feature.gff"
printf "chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=transcript1\n" >> "$TMPDIR/feature.gff"
printf "chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tID=cds1;Parent=transcript1\n" >> "$TMPDIR/feature.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "$TMPDIR/feature.gff"
printf "chr2\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "$TMPDIR/feature.gff"
printf "chr3\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "$TMPDIR/feature.gff"

# Create expected output files
printf "chr1\t100\t250\nchr1\t300\t400\n" > "$TMPDIR/expected.bed"
printf "chr1\t100\t250\nchr1\t250\t500\nchr1\t501\t1000\n" > "$TMPDIR/expected_strand.bed"
printf "chr1\t100\t250\nchr1\t501\t1000\n" > "$TMPDIR/expected_specific_strand.bed"
printf "chr1\t128\t228\nchr1\t428\t528\n" > "$TMPDIR/expected_bam.bed"
printf "chr1\t100\t400\n" > "$TMPDIR/expected_distance.bed"
printf "chr1\t100\t500\t2\t1\t3\nchr1\t501\t1000\t4\t4\t4\n" > "$TMPDIR/expected_operation.bed"
printf "chr1\t100\t500\ta1|a2|a3\nchr1\t501\t1000\ta4\n" > "$TMPDIR/expected_delim.bed"
printf "chr1\t100\t500\t2.567\nchr1\t501\t1000\t4\n" > "$TMPDIR/expected_precision.bed"
printf "##gff-version 3\nchr1\t999\t2000\nchr2\t1499\t1700\nchr3\t999\t2000\n" > "$TMPDIR/expected_header.bed"

# Test 1: Default sort on BED file
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bedtools_merge on BED file"
"$meta_executable" \
  --input "../featureA.bed" \
  --output "output.bed"

# # checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected.bed"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2: strand option
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bedtools_merge on BED file with strand option"
"$meta_executable" \
  --input "../featureB.bed" \
  --output "output.bed" \
  --strand

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_strand.bed"
echo "- test2 succeeded -"

popd > /dev/null

# Test 3: specific strand option
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null

echo "> Run bedtools_merge on BED file with specific strand option"
"$meta_executable" \
  --input "../featureB.bed" \
  --output "output.bed" \
  --specific_strand "+" 

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_specific_strand.bed"
echo "- test3 succeeded -"

popd > /dev/null

# Test 4: BED option
mkdir "$TMPDIR/test4" && pushd "$TMPDIR/test4" > /dev/null

echo "> Run bedtools_merge on BAM file with BED option"
"$meta_executable" \
  --input "$test_data/feature.bam" \
  --output "output.bed" \
  --bed

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_bam.bed"
echo "- test4 succeeded -"

popd > /dev/null

# Test 5: distance option
mkdir "$TMPDIR/test5" && pushd "$TMPDIR/test5" > /dev/null

echo "> Run bedtools_merge on BED file with distance option"
"$meta_executable" \
  --input "../featureA.bed" \
  --output "output.bed" \
  --distance -5 

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected.bed"
echo "- test5 succeeded -"

popd > /dev/null

# Test 6: columns option & operation option
mkdir "$TMPDIR/test6" && pushd "$TMPDIR/test6" > /dev/null

echo "> Run bedtools_merge on BED file with columns & operation options"
"$meta_executable" \
  --input "../featureB.bed" \
  --output "output.bed" \
  --columns 5 \
  --operation "mean,min,max"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_operation.bed"
echo "- test6 succeeded -"

popd > /dev/null

# Test 7: delimeter option
mkdir "$TMPDIR/test7" && pushd "$TMPDIR/test7" > /dev/null

echo "> Run bedtools_merge on BED file with delimeter option"
"$meta_executable" \
  --input "../featureB.bed" \
  --output "output.bed" \
  --columns 4 \
  --operation "collapse" \
  --delimiter "|"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_delim.bed"
echo "- test7 succeeded -"

popd > /dev/null

# Test 8: precision option
mkdir "$TMPDIR/test8" && pushd "$TMPDIR/test8" > /dev/null

echo "> Run bedtools_merge on BED file with precision option"
"$meta_executable" \
  --input "../feature_precision.bed" \
  --output "output.bed" \
  --columns 5 \
  --operation "mean" \
  --precision 4

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_precision.bed"
echo "- test8 succeeded -"

popd > /dev/null

# Test 9: header option
mkdir "$TMPDIR/test9" && pushd "$TMPDIR/test9" > /dev/null

echo "> Run bedtools_merge on GFF file with header option"
"$meta_executable" \
  --input "../feature.gff" \
  --output "output.gff" \
  --header

# checks
assert_file_exists "output.gff"
assert_file_not_empty "output.gff"
assert_identical_content "output.gff" "../expected_header.bed"
echo "- test9 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0
