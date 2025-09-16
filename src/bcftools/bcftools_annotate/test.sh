#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Source test helpers
source "$meta_resources_dir/test_helpers.sh"

# Create directories for tests
log "Creating Test Data..."
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

bgzip -c $TMPDIR/example.vcf > $TMPDIR/example.vcf.gz
tabix -p vcf $TMPDIR/example.vcf.gz

cat <<EOF > "$TMPDIR/annots.tsv"
1	752567	752567	FooValue1	12345
1	752722	752722	FooValue2	67890
EOF

cat <<EOF > "$TMPDIR/rename.tsv"
INFO/.	Luigi
EOF

bgzip $TMPDIR/annots.tsv
tabix -s1 -b2 -e3 $TMPDIR/annots.tsv.gz

cat <<EOF > "$TMPDIR/header.hdr"
##FORMAT=<ID=FOO,Number=1,Type=String,Description="Some description">
##INFO=<ID=BAR,Number=1,Type=Integer,Description="Some description">
EOF

cat <<EOF > "$TMPDIR/rename_chrm.tsv"
1	chr1
2	chr2
EOF

# Test 1: Remove ID annotations
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

log "Test 1: Remove ID annotations"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --remove "ID" \

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "1	752567	.	A	C" "VCF with removed ID"
log "✓ Test 1 passed"

popd > /dev/null

# Test 2: Annotate with -a, -c and -h options
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

log "Test 2: Annotate with -a, -c and -h options"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --annotations "../annots.tsv.gz" \
  --header_lines "../header.hdr" \
  --columns "CHROM,FROM,TO,FMT/FOO,BAR" \
  --mark_sites "BAR" \

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" $(echo -e "1\t752567\tllama\tA\tC\t.\t.\tBAR=12345\tFOO\tFooValue1") "annotated VCF content"
log "✓ Test 2 passed"

popd > /dev/null

# Test 3: Set ID option
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null

log "Test 3: Set ID option"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --set_id "+'%CHROM\_%POS\_%REF\_%FIRST_ALT'" \

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "'1_752722_G_A'" "VCF with set ID"
log "✓ Test 3 passed"

popd > /dev/null

# Test 4: Rename annotations
mkdir "$TMPDIR/test4" && pushd "$TMPDIR/test4" > /dev/null

log "Test 4: Rename annotations"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --rename_annotations "../rename.tsv"

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "##bcftools_annotateCommand=annotate --rename-annots ../rename.tsv -o annotated.vcf" "VCF with command line"
log "✓ Test 4 passed"

popd > /dev/null

# Test 5: Rename chromosomes
mkdir "$TMPDIR/test5" && pushd "$TMPDIR/test5" > /dev/null

log "Test 5: Rename chromosomes"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --rename_chromosomes "../rename_chrm.tsv"

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "chr1" "VCF with renamed chromosomes"
log "✓ Test 5 passed"

popd > /dev/null

# Test 6: Sample option
mkdir "$TMPDIR/test6" && pushd "$TMPDIR/test6" > /dev/null

log "Test 6: Sample selection"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --samples "SAMPLE1"

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "##bcftools_annotateCommand=annotate -s SAMPLE1 -o annotated.vcf ../example.vcf" "VCF with sample selection"
log "✓ Test 6 passed"

popd > /dev/null

# Test 7: Single overlaps
mkdir "$TMPDIR/test7" && pushd "$TMPDIR/test7" > /dev/null

log "Test 7: Single overlaps and keep sites"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --single_overlaps \
  --keep_sites \

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "annotate -k --single-overlaps -o annotated.vcf ../example.vcf" "VCF with single overlaps option"
log "✓ Test 7 passed"

popd > /dev/null

# Test 8: Min overlap
mkdir "$TMPDIR/test8" && pushd "$TMPDIR/test8" > /dev/null

log "Test 8: Minimum overlap option"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --annotations "../annots.tsv.gz" \
  --columns "CHROM,FROM,TO,FMT/FOO,BAR" \
  --header_lines "../header.hdr" \
  --min_overlap "1"

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "annotate -a ../annots.tsv.gz -c CHROM,FROM,TO,FMT/FOO,BAR -h ../header.hdr --min-overlap 1 -o annotated.vcf ../example.vcf" "VCF with min overlap"
log "✓ Test 8 passed"

popd > /dev/null

# Test 9: Regions
mkdir "$TMPDIR/test9" && pushd "$TMPDIR/test9" > /dev/null

log "Test 9: Region filtering"
"$meta_executable" \
  --input "../example.vcf.gz" \
  --output "annotated.vcf" \
  --regions "1:752567-752722"

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "annotate -r 1:752567-752722 -o annotated.vcf ../example.vcf.gz" "VCF with region filtering"
log "✓ Test 9 passed"

popd > /dev/null

# Test 10: Pair logic
mkdir "$TMPDIR/test10" && pushd "$TMPDIR/test10" > /dev/null

log "Test 10: Pair logic option"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --pair_logic "all"

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "annotate --pair-logic all -o annotated.vcf ../example.vcf" "VCF with pair logic"
log "✓ Test 10 passed"

popd > /dev/null

# Test 11: Regions overlap
mkdir "$TMPDIR/test11" && pushd "$TMPDIR/test11" > /dev/null

log "Test 11: Regions overlap option"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --regions_overlap "1"

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "annotate --regions-overlap 1 -o annotated.vcf ../example.vcf" "VCF with regions overlap"
log "✓ Test 11 passed"

popd > /dev/null

# Test 12: Include filter
mkdir "$TMPDIR/test12" && pushd "$TMPDIR/test12" > /dev/null

log "Test 12: Include filter expression"
"$meta_executable" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --include "FILTER='PASS'" \

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "annotate -i FILTER='PASS' -o annotated.vcf ../example.vcf" "VCF with include filter"
log "✓ Test 12 passed"

popd > /dev/null

# Test 13: Exclude filter with merge logic
mkdir "$TMPDIR/test13" && pushd "$TMPDIR/test13" > /dev/null

log "Test 13: Exclude filter with merge logic"
"$meta_executable" \
  --annotations "../annots.tsv.gz" \
  --input "../example.vcf" \
  --output "annotated.vcf" \
  --exclude "FILTER='PASS'" \
  --header_lines "../header.hdr" \
  --columns "CHROM,FROM,TO,FMT/FOO,BAR" \
  --merge_logic "FOO:first" \

# checks
check_file_exists "annotated.vcf" "output VCF file"
check_file_not_empty "annotated.vcf" "output VCF file"
check_file_contains "annotated.vcf" "annotate -a ../annots.tsv.gz -c CHROM,FROM,TO,FMT/FOO,BAR -e FILTER='PASS' -h ../header.hdr -l FOO:first -o annotated.vcf ../example.vcf" "VCF with exclude filter and merge logic"
log "✓ Test 13 passed"

popd > /dev/null

log "All tests completed successfully!"
exit 0

