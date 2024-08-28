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
1	752567	llama	A	C	.	.	.	.	.
1	752722	.	G	A	.	.	.	.	.
EOF

cat <<EOF > "$TMPDIR/annots.tsv"
1	752567	752567	FooValue1	12345
1	752722	752722	FooValue2	67890
EOF

bgzip $TMPDIR/annots.tsv
tabix -s1 -b2 -e3 $TMPDIR/annots.tsv.gz

cat <<EOF > "$TMPDIR/header.hdr"
##FORMAT=<ID=FOO,Number=1,Type=String,Description="Some description">
##INFO=<ID=BAR,Number=1,Type=Integer,Description="Some description">
EOF

# Test 1: Remove ID annotations
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bcftools_annotate remove annotations"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --remove "ID"

# checks
assert_file_exists "annotated.vcf"
assert_file_not_empty "annotated.vcf"
assert_file_contains "annotated.vcf" "1	752567	.	A	C"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2: Annotate with -a, -c and -h options
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bcftools_annotate with -a, -c and -h options"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --annotations "../annots.tsv.gz" \
  --header_lines "../header.hdr" \
  --columns "CHROM,FROM,TO,FMT/FOO,BAR"

# checks
assert_file_exists "annotated.vcf"
assert_file_not_empty "annotated.vcf"
assert_file_contains "annotated.vcf" $(echo -e "1\t752567\tllama\tA\tC\t.\t.\tBAR=12345\tFOO\tFooValue1")
echo "- test2 succeeded -"

popd > /dev/null

# Test 3: 
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null

echo "> Run bcftools_annotate with --set_id option"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --set_id "+'%CHROM\_%POS\_%REF\_%FIRST_ALT'" \

# checks
assert_file_exists "annotated.vcf"
assert_file_not_empty "annotated.vcf"
assert_file_contains "annotated.vcf" "'1_752722_G_A'"
echo "- test3 succeeded -"

popd > /dev/null

# echo
# cat "annotated.vcf"
# echo
# exit 0

echo "---- All tests succeeded! ----"
exit 0