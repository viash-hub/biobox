#!/bin/bash

set -e

#############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_doesnt_exist() {
  [ ! -f "$1" ] || { echo "File '$1' exists but shouldn't" && exit 1; }
}
assert_file_empty() {
  [ ! -s "$1" ] || { echo "File '$1' is not empty but should be" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_not_contains() {
  grep -q "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}
#############################################

in_fa="$meta_resources_dir/test_data/reference_small.fa"
in_gtf="$meta_resources_dir/test_data/reference_small.gtf"

echo "#############################################"
echo "> Simple run"

mkdir simple_run
cd simple_run

out_tar="myreference.tar.gz"

echo "> Running $meta_name."
$meta_executable \
  --genome_fasta "$in_fa" \
  --gtf "$in_gtf" \
  --reference_archive "$out_tar" \
  --extra_star_params "--genomeSAindexNbases 6" \
  ---cpus 2

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

assert_file_exists "$out_tar"
assert_file_not_empty "$out_tar"

echo ">> Checking whether output contains the expected files"
tar -xvf "$out_tar" > /dev/null
assert_file_exists "BD_Rhapsody_Reference_Files/star_index/genomeParameters.txt"
assert_file_exists "BD_Rhapsody_Reference_Files/bwa-mem2_index/reference_small.ann"
assert_file_exists "BD_Rhapsody_Reference_Files/reference_small-processed.gtf"
assert_file_exists "BD_Rhapsody_Reference_Files/mitochondrial_contigs.txt"
assert_file_contains "BD_Rhapsody_Reference_Files/reference_small-processed.gtf" "chr1.*HAVANA.*ENSG00000243485"
assert_file_contains "BD_Rhapsody_Reference_Files/mitochondrial_contigs.txt" 'chrMT'

cd ..

echo "#############################################"

echo "> Tests succeeded!"