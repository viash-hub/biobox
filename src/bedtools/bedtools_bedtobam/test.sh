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

# Expected output sam files for each test
cat <<EOF > "$TMPDIR/expected.sam"
@HD	VN:1.0	SO:unsorted
@PG	ID:BEDTools_bedToBam	VN:Vv2.30.0
@PG	ID:samtools	PN:samtools	PP:BEDTools_bedToBam	VN:1.16.1	CL:samtools view -h output.bam
@SQ	SN:chr1	AS:../genome.txt	LN:248956422
@SQ	SN:chr3	AS:../genome.txt	LN:242193529
@SQ	SN:chr2	AS:../genome.txt	LN:198295559
my_read/1	0	chr1	129	255	100M	*	0	0	*	*
my_read/2	16	chr1	429	255	100M	*	0	0	*	*
EOF
cat <<EOF > "$TMPDIR/expected12.sam"
@HD	VN:1.0	SO:unsorted
@PG	ID:BEDTools_bedToBam	VN:Vv2.30.0
@PG	ID:samtools	PN:samtools	PP:BEDTools_bedToBam	VN:1.16.1	CL:samtools view -h output.bam
@SQ	SN:chr1	AS:../genome.txt	LN:248956422
@SQ	SN:chr3	AS:../genome.txt	LN:242193529
@SQ	SN:chr2	AS:../genome.txt	LN:198295559
my_read/1	0	chr1	129	255	100M	*	0	0	*	*
my_read/2	16	chr1	429	255	100M	*	0	0	*	*
EOF
cat <<EOF > "$TMPDIR/expected_mapquality.sam"
@HD	VN:1.0	SO:unsorted
@PG	ID:BEDTools_bedToBam	VN:Vv2.30.0
@PG	ID:samtools	PN:samtools	PP:BEDTools_bedToBam	VN:1.16.1	CL:samtools view -h output.bam
@SQ	SN:chr1	AS:../genome.txt	LN:248956422
@SQ	SN:chr3	AS:../genome.txt	LN:242193529
@SQ	SN:chr2	AS:../genome.txt	LN:198295559
my_read/1	0	chr1	129	10	100M	*	0	0	*	*
my_read/2	16	chr1	429	10	100M	*	0	0	*	*
EOF
cat <<EOF > "$TMPDIR/expected_gff.sam"
@HD	VN:1.0	SO:unsorted
@PG	ID:BEDTools_bedToBam	VN:Vv2.30.0
@PG	ID:samtools	PN:samtools	PP:BEDTools_bedToBam	VN:1.16.1	CL:samtools view -h output.bam
@SQ	SN:chr1	AS:../genome.txt	LN:248956422
@SQ	SN:chr3	AS:../genome.txt	LN:242193529
@SQ	SN:chr2	AS:../genome.txt	LN:198295559
gene	0	chr1	1000	255	1001M	*	0	0	*	*
mRNA	0	chr3	1000	255	1001M	*	0	0	*	*
exon	0	chr1	1000	255	201M	*	0	0	*	*
exon	0	chr2	1500	255	201M	*	0	0	*	*
CDS	0	chr1	1000	255	201M	*	0	0	*	*
CDS	0	chr1	1500	255	201M	*	0	0	*	*
EOF

# Test 1: Default conversion BED to BAM
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bedtools_bedtobam on BED file"
"$meta_executable" \
  --input "../example.bed" \
  --genome "../genome.txt" \
  --output "output.bam"

samtools view -h output.bam > output.sam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
assert_identical_content "output.sam" "../expected.sam"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2: BED12 file
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bedtools_bedtobam on BED12 file"
"$meta_executable" \
  --input "../example.bed12" \
  --genome "../genome.txt" \
  --output "output.bam" \
  --bed12 \

samtools view -h output.bam > output.sam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
assert_identical_content "output.sam" "../expected12.sam"
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
# Cannot assert_identical_content because umcompress option does not work on this version of bedtools.

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

samtools view -h output.bam > output.sam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
assert_identical_content "output.sam" "../expected_mapquality.sam"
echo "- test4 succeeded -"

popd > /dev/null

# Test 5: gff to bam conversion
mkdir "$TMPDIR/test5" && pushd "$TMPDIR/test5" > /dev/null

echo "> Run bedtools_bedtobam on GFF file"
"$meta_executable" \
  --input "../example.gff" \
  --genome "../genome.txt" \
  --output "output.bam"

samtools view -h output.bam > output.sam

# checks
assert_file_exists "output.bam"
assert_file_not_empty "output.bam"
assert_identical_content "output.sam" "../expected_gff.sam"
echo "- test5 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0
