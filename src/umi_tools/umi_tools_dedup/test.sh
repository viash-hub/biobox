#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/out"

mkdir -p "$out_dir"

############################################################################################

echo ">>> Test 1: Basic usage of $meta_functionality_name"

"$meta_executable" \
  --paired \
  --input "$test_dir/sample.bam" \
  --bai "$test_dir/sample.bam.bai" \
  --output "$out_dir/deduped.bam"

echo ">>> Checking whether output exists"
[ ! -f "$out_dir/deduped.bam" ] && echo "File 'deduped.bam' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$out_dir/deduped.bam" ] && echo "File 'deduped.bam' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$out_dir/deduped.bam" "$test_dir/deduped.bam" || \
    (echo "Output file deduped.bam does not match expected output" && exit 1)


############################################################################################

rm -rf "$out_dir"

echo "All tests succeeded!"
exit 0