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
TMPDIR=$(mktemp -d "$meta_temp_dir/XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Create and populate input files
printf "chr1\t248956422\nchr3\t242193529\nchr2\t198295559\n" > "$TMPDIR/genome.txt"
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t37\t+\nchr2:172936693-172938111\t428\t528\tmy_read/2\t37\t-\n" > "$TMPDIR/example.bed"
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\t128\t228\t255,0,0\t1\t100\t0\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\t428\t528\t255,0,0\t1\t100\t0\n" > "$TMPDIR/example.bed12"
# Create and populate example.gff file
printf "##gff-version 3\n" > "$TMPDIR/example.gff"
printf "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "$TMPDIR/example.gff"
printf "chr3\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "$TMPDIR/example.gff"
printf "chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=transcript1\n" >> "$TMPDIR/example.gff"
printf "chr2\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "$TMPDIR/example.gff"
printf "chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tID=cds1;Parent=transcript1\n" >> "$TMPDIR/example.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "$TMPDIR/example.gff"

# Test 1: Default conversion BED to BAM
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bedtools_bedtobam on BED file"
"$meta_executable" \
  --input "../example.bed" \
  --genome "../genome.txt" \
  --output "output.bam"

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
echo "- test1 succeeded -"

popd > /dev/null

# I could not assert_identical_content (bam file is compressed and -ubam option does not seem to work on this version).
# could use samtools view to check the content of the bam file, would need to add samtools as a dependency.
# TODO: check if older versions works with -ubam option!!!
# Version 2.27.1 (bedtools) -ubam outputs better than version 2.30.0 (bedtools).
# However, its output is still somewhat compressed.

# Test 2: BED12 file
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bedtools_bedtobam on BED12 file"
"$meta_executable" \
  --input "../example.bed12" \
  --genome "../genome.txt" \
  --output "output.bam" \
  --bed12 \
#  --uncompress_bam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
# assert_identical_content "output.bam" "../test_data/expected12.sam"
echo "- test2 succeeded -"

popd > /dev/null

# Test 3: Uncompressed BAM file
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null

echo "> Run bedtools_bedtobam on BED file with uncompressed BAM output"
"$meta_executable" \
  --input "../example.bed" \
  --genome "../genome.txt" \
  --output "output.bam" \
  --uncompress_bam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
# assert_identical_content "output.bam" "../test_data/expected.sam"

echo "- test3 succeeded -"

popd > /dev/null

# Test 4: Map quality
mkdir "$TMPDIR/test4" && pushd "$TMPDIR/test4" > /dev/null

echo "> Run bedtools_bedtobam on BED file with map quality"
"$meta_executable" \
  --input "../example.bed" \
  --genome "../genome.txt" \
  --output "output.bam" \
  --map_quality 10 \
#  --uncompress_bam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
#assert_identical_content "output.bam" "../test_data/expected_mapquality.sam"
echo "- test4 succeeded -"

popd > /dev/null

# Test 5: gff to bam conversion
mkdir "$TMPDIR/test5" && pushd "$TMPDIR/test5" > /dev/null

echo "> Run bedtools_bedtobam on GFF file"
"$meta_executable" \
  --input "../example.gff" \
  --genome "../genome.txt" \
  --output "output.bam"
#  --uncompress_bam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
# assert_identical_content "output.bam" "../test_data/expected_gff.sam"
echo "- test5 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0
