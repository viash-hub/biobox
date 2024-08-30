#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

#test_data="$meta_resources_dir/test_data"

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

# Create test data
cat <<EOF > "$TMPDIR/example.vcf"
##fileformat=VCFv4.1
##contig=<ID=1,length=249250621,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	752567	llama	G	C,A	.	.	.	.	1/2
1	752752	.	G	A,AAA	.	.	.	.	./.
EOF

bgzip -c $TMPDIR/example.vcf > $TMPDIR/example.vcf.gz
tabix -p vcf $TMPDIR/example.vcf.gz

cat <<EOF > "$TMPDIR/example_2.vcf"
##fileformat=VCFv4.1
##contig=<ID=1,length=249250621,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	752569	cat	G	C,A	.	.	.	.	1/2
1	752739	.	G	A,AAA	.	.	.	.	./.
EOF

bgzip -c $TMPDIR/example_2.vcf > $TMPDIR/example_2.vcf.gz
tabix -p vcf $TMPDIR/example_2.vcf.gz

cat <<EOF > "$TMPDIR/file_list.txt"
$TMPDIR/example.vcf.gz
$TMPDIR/example_2.vcf.gz
EOF

# Test 1: Default test
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bcftools_concat default test"
"$meta_executable" \
  --input "../example.vcf" \
  --input "../example_2.vcf" \
  --output "concatenated.vcf" \

#  &> /dev/null

# checks
assert_file_exists "concatenated.vcf"
assert_file_not_empty "concatenated.vcf"
assert_file_contains "concatenated.vcf" "concat -o concatenated.vcf ../example.vcf ../example_2.vcf"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2: Allow overlaps, compact PS and remove duplicates
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bcftools_concat test with allow overlaps, and remove duplicates"
"$meta_executable" \
  --input "../example.vcf.gz" \
  --input "../example_2.vcf.gz" \
  --output "concatenated.vcf" \
  --allow_overlaps \
  --remove_duplicates 'none' \
#  &> /dev/null

# checks
assert_file_exists "concatenated.vcf"
assert_file_not_empty "concatenated.vcf"
assert_file_contains "concatenated.vcf" "concat -a -d none -o concatenated.vcf ../example.vcf.gz ../example_2.vcf.gz"  
echo "- test2 succeeded -"

popd > /dev/null


# Test 3: Ligate, ligate force and ligate warn
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null

echo "> Run bcftools_concat test with ligate, ligate force and ligate warn"
"$meta_executable" \
  --input "../example.vcf.gz" \
  --input "../example_2.vcf.gz" \
  --output "concatenated.vcf" \
  --ligate \
  --compact_PS \
#  &> /dev/null


# checks
assert_file_exists "concatenated.vcf"
assert_file_not_empty "concatenated.vcf"
assert_file_contains "concatenated.vcf" "concat -c -l -o concatenated.vcf ../example.vcf.gz ../example_2.vcf.gz"
echo "- test3 succeeded -"

popd > /dev/null

# Test 4: file list with ligate force and ligate warn
mkdir "$TMPDIR/test4" && pushd "$TMPDIR/test4" > /dev/null

echo "> Run bcftools_concat test with file list, ligate force and ligate warn"
"$meta_executable" \
  --file_list "../file_list.txt" \
  --output "concatenated.vcf" \
  --ligate_force \
#  &> /dev/null

# checks
assert_file_exists "concatenated.vcf"
assert_file_not_empty "concatenated.vcf"
assert_file_contains "concatenated.vcf" "concat --ligate-force -o concatenated.vcf -f ../file_list.txt"
echo "- test4 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0

echo
cat concatenated.vcf
echo

