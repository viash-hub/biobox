#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/bedtools/bedtools_intersect/bedtools_intersect"
meta_resources_dir="src/bedtools/bedtools_intersect"
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

# Create and populate featuresA.bed
printf "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "test_data/featuresA.bed"

# Create and populate featuresB.bed
printf "chr1\t180\t280\nchr1\t290\t390\nchr1\t500\t600\n" > "test_data/featuresB.bed"

# Create and populate featuresC.bed
printf "chr1\t120\t220\nchr1\t250\t350\nchr1\t500\t580\n" > "test_data/featuresC.bed"

# Create and populate examples gff files
# example1.gff
printf "##gff-version 3\n" > "test_data/example1.gff"
printf "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/example1.gff"
printf "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/example1.gff"
printf "chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=transcript1\n" >> "test_data/example1.gff"
printf "chr1\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/example1.gff"
printf "chr1\t.\tCDS\t1000\t1200\t.\t+\t0\tID=cds1;Parent=transcript1\n" >> "test_data/example1.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/example1.gff"
# example2.gff
printf "##gff-version 3\n" > "test_data/example2.gff"
printf "chr1\t.\tgene\t1200\t1800\t.\t-\t.\tID=gene2;Name=Gene2\n" >> "test_data/example2.gff"
printf "chr1\t.\tmRNA\t1400\t2000\t.\t-\t.\tID=transcript2;Parent=gene2\n" >> "test_data/example2.gff"
printf "chr1\t.\texon\t1400\t2000\t.\t-\t.\tID=exon3;Parent=transcript2\n" >> "test_data/example2.gff"
printf "chr1\t.\texon\t1600\t2000\t.\t-\t.\tID=exon4;Parent=transcript2\n" >> "test_data/example2.gff"
printf "chr1\t.\tCDS\t3000\t3200\t.\t-\t1\tID=cds3;Parent=transcript2\n" >> "test_data/example2.gff"
printf "chr1\t.\tCDS\t3500\t3700\t.\t-\t0\tID=cds4;Parent=transcript2\n" >> "test_data/example2.gff"

# Create and populate expected output files for different tests
printf "chr1\t180\t200\nchr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_default.bed"
printf "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "test_data/expected_wa.bed"
printf "chr1\t180\t200\tchr1\t180\t280\nchr1\t180\t250\tchr1\t180\t280\nchr1\t300\t390\tchr1\t290\t390\n" > "test_data/expected_wb.bed"
printf "chr1\t100\t200\tchr1\t180\t280\nchr1\t150\t250\tchr1\t180\t280\nchr1\t300\t400\tchr1\t290\t390\n" > "test_data/expected_loj.bed"
printf "chr1\t100\t200\tchr1\t180\t280\t20\nchr1\t150\t250\tchr1\t180\t280\t70\nchr1\t300\t400\tchr1\t290\t390\t90\n" > "test_data/expected_wo.bed"
printf "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "test_data/expected_u.bed"
printf "chr1\t100\t200\t1\nchr1\t150\t250\t1\nchr1\t300\t400\t1\n" > "test_data/expected_c.bed"
printf "chr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_f50.bed"
printf "chr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_f30.bed"
printf "chr1\t180\t200\nchr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_f10.bed"
printf "chr1\t180\t200\nchr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_r.bed"
printf "chr1\t180\t200\nchr1\t120\t200\nchr1\t180\t250\nchr1\t150\t220\nchr1\t300\t390\nchr1\t300\t350\n" > "test_data/expected_multiple.bed"
# expected gff file
printf "chr1\t.\tgene\t1200\t1800\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/expected.gff"
printf "chr1\t.\tgene\t1400\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/expected.gff"
printf "chr1\t.\tgene\t1400\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/expected.gff"
printf "chr1\t.\tgene\t1600\t2000\t.\t+\t.\tID=gene1;Name=Gene1\n" >> "test_data/expected.gff"
printf "chr1\t.\tmRNA\t1200\t1800\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/expected.gff"
printf "chr1\t.\tmRNA\t1400\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/expected.gff"
printf "chr1\t.\tmRNA\t1400\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/expected.gff"
printf "chr1\t.\tmRNA\t1600\t2000\t.\t+\t.\tID=transcript1;Parent=gene1\n" >> "test_data/expected.gff"
printf "chr1\t.\texon\t1200\t1200\t.\t+\t.\tID=exon1;Parent=transcript1\n" >> "test_data/expected.gff"
printf "chr1\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/expected.gff"
printf "chr1\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/expected.gff"
printf "chr1\t.\texon\t1500\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/expected.gff"
printf "chr1\t.\texon\t1600\t1700\t.\t+\t.\tID=exon2;Parent=transcript1\n" >> "test_data/expected.gff"
printf "chr1\t.\tCDS\t1200\t1200\t.\t+\t0\tID=cds1;Parent=transcript1\n" >> "test_data/expected.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/expected.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/expected.gff"
printf "chr1\t.\tCDS\t1500\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/expected.gff"
printf "chr1\t.\tCDS\t1600\t1700\t.\t+\t2\tID=cds2;Parent=transcript1\n" >> "test_data/expected.gff"

# Test 1: Default intersect
mkdir test1
cd test1

echo "> Run bedtools_intersect on BED files with default intersect"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_default.bed"
echo "- test1 succeeded -"

cd ..

# Test 2: Write A option
mkdir test2
cd test2

echo "> Run bedtools_intersect on BED files with -wa option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --write_a

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_wa.bed"
echo "- test2 succeeded -"

cd ..

# Test 3: -wb option
mkdir test3
cd test3

echo "> Run bedtools_intersect on BED files with -wb option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --write_b 

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_wb.bed"
echo "- test3 succeeded -"

cd ..

# Test 4: -loj option
mkdir test4
cd test4

echo "> Run bedtools_intersect on BED files with -loj option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --left_outer_join

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_loj.bed"
echo "- test4 succeeded -"

cd ..

# Test 5: -wo option
mkdir test5
cd test5

echo "> Run bedtools_intersect on BED files with -wo option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --write_overlap 
  

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_wo.bed"
echo "- test5 succeeded -"

cd ..

# Test 6: -u option
mkdir test6
cd test6

echo "> Run bedtools_intersect on BED files with -u option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --report_A_if_no_overlap true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_u.bed"
echo "- test6 succeeded -"

cd ..

# Test 7: -c option
mkdir test7
cd test7

echo "> Run bedtools_intersect on BED files with -c option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --number_of_overlaps_A true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_c.bed"
echo "- test7 succeeded -"

cd ..

# Test 8: -f 0.50 option
mkdir test8
cd test8

echo "> Run bedtools_intersect on BED files with -f 0.50 option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --min_overlap_A 0.50

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_f50.bed"
echo "- test8 succeeded -"

cd ..

# Test 9: -f 0.30 option
mkdir test9
cd test9

echo "> Run bedtools_intersect on BED files with -f 0.30 option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --min_overlap_A 0.30

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_f30.bed"
echo "- test9 succeeded -"

cd ..

# Test 10: -f 0.10 option
mkdir test10
cd test10

echo "> Run bedtools_intersect on BED files with -f 0.10 option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --min_overlap_A 0.10

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_f10.bed"
echo "- test10 succeeded -"

cd ..

# Test 11: -r option
mkdir test11
cd test11

echo "> Run bedtools_intersect on BED files with -r option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --reciprocal_overlap true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_r.bed"
echo "- test11 succeeded -"

cd ..


# Test 12: Multiple files
mkdir test12
cd test12

echo "> Run bedtools_intersect on multiple BED files"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --input_b "../test_data/featuresC.bed" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_multiple.bed"
echo "- test12 succeeded -"

cd ..

# Test 13: VCF file format
mkdir test13
cd test13

echo "> Run bedtools_intersect on GFF files"
"$meta_executable" \
  --input_a "../test_data/example1.gff" \
  --input_b "../test_data/example2.gff" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected.gff"
echo "- test13 succeeded -"

cd ..

echo "---- All tests succeeded! ----"
exit 0
