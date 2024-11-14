#!/bin/bash

set -eou pipefail

# Helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}


# create temporary directory and clean up on exit
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXX")
function clean_up {
 [[ -d "$TMPDIR" ]] && rm -rf "$TMPDIR"
}
trap clean_up EXIT

output_test1="$TMPDIR/output1"
mkdir "$output_test1"
sample_dir_test_1="$output_test1/fastq"
mkdir "$sample_dir_test_1"

"$meta_executable" \
    --fastqs "$meta_resources_dir/test_data/fastq" \
    --sample_metadata "$meta_resources_dir/test_data/samplesheet.csv" \
    --sample_fastq "$sample_dir_test_1" \
    --metrics "$output_test1/metrics.tsv" \
    --most_frequent_unmatched "$output_test1/most_frequent_unmatched.tsv" \
    --sample_barcode_hop_metrics "$output_test1/sample_barcode_hop_metrics.tsv" \
    --per_sample_metrics "$output_test1/per_sample_metrics.tsv" \
    --per_project_metrics "$output_test1/per_project_metrics.tsv" \
    ---cpus 1
    
# Check for correct number of output FASTQ files
readarray -d '' output_fastq < <(find "$sample_dir_test_1" -name "*.fastq.gz" -print0)
if (( ${#output_fastq[@]} != "196" )); then
    echo "Wrong number of output fastq files found."
    exit 1
fi

# Check if fastq files are not empty
for fastq in ${output_fastq[@]}; do
   assert_file_not_empty "$fastq"
done

# Checking if requested output files exist
assert_file_exists "$output_test1/metrics.tsv"
assert_file_exists "$output_test1/most_frequent_unmatched.tsv"
assert_file_exists "$output_test1/sample_barcode_hop_metrics.tsv"
assert_file_exists "$output_test1/per_sample_metrics.tsv"
assert_file_exists "$output_test1/per_project_metrics.tsv"

# Checking output file contents
diff -q "$meta_resources_dir/test_data/expected/metrics.tsv" "$output_test1/metrics.tsv" || \
    (echo "Incorrect metrics.tsv output!" && exit 1)

diff -q "$meta_resources_dir/test_data/expected/per_project_metrics.tsv" "$output_test1/per_project_metrics.tsv" || \
  (echo "Incorrect per_project_metrics.tsv output!" && diff exit 1)

diff -q "$meta_resources_dir/test_data/expected/per_sample_metrics.tsv" "$output_test1/per_sample_metrics.tsv" || \
  (echo "Incorrect per_sample_metrics.tsv output!" && exit 1) 