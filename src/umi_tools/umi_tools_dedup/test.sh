#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/out"

mkdir -p "$out_dir"

############################################################################################

echo ">>> Test 1: Basic usage of $meta_functionality_name with statistics output"

"$meta_executable" \
  --paired \
  --input "$test_dir/sample.bam" \
  --bai "$test_dir/sample.bam.bai" \
  --output "$out_dir/deduped.sam" \
  --out_sam \
  --output_stats "$out_dir/dedup" \
  --random_seed 1

echo ">>> Checking whether output exists"
[ ! -f "$out_dir/deduped.sam" ] && echo "File 'deduped.sam' does not exist!" && exit 1
[ ! -f "$out_dir/dedup_edit_distance.tsv" ] && echo "File 'dedup_edit_distance.tsv' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$out_dir/deduped.sam" ] && echo "File 'deduped.sam' is empty!" && exit 1
[ ! -s "$out_dir/dedup_edit_distance.tsv" ] && echo "File 'dedup_edit_distance.tsv' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$out_dir/deduped.sam" "$test_dir/deduped.sam" || \
    (echo "Output file deduped.sam does not match expected output" && exit 1)
diff "$out_dir/dedup_edit_distance.tsv" "$test_dir/dedup_edit_distance.tsv" || \
    (echo "Output file dedup_edit_distance.tsv does not match expected output" && exit 1)

############################################################################################

echo ">>> Test 2: $meta_functionality_name with random subset selection"

"$meta_executable" \
  --paired \
  --input "$test_dir/sample.bam" \
  --bai "$test_dir/sample.bam.bai" \
  --output "$out_dir/deduped_fraction.sam" \
  --out_sam \
  --subset 0.5 \
  --random_seed 1


echo ">>> Checking whether output exists"
[ ! -f "$out_dir/deduped_fraction.sam" ] && echo "File 'deduped_fraction.sam' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$out_dir/deduped_fraction.sam" ] && echo "File 'deduped_fraction.sam' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$out_dir/deduped_fraction.sam" "$test_dir/deduped_fraction.sam" || \
    (echo "Output file deduped_fraction.sam does not match expected output" && exit 1)

############################################################################################

echo ">>> Test 3: $meta_functionality_name with --method unique"

"$meta_executable" \
  --paired \
  --input "$test_dir/sample.bam" \
  --bai "$test_dir/sample.bam.bai" \
  --output "$out_dir/deduped_unique.sam" \
  --out_sam \
  --method "unique" \
  --random_seed 1

echo ">>> Checking whether output exists"
[ ! -f "$out_dir/deduped_unique.sam" ] && echo "File 'deduped_unique.sam' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$out_dir/deduped_unique.sam" ] && echo "File 'deduped_unique.sam' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$out_dir/deduped_unique.sam" "$test_dir/deduped_unique.sam" || \
    (echo "Output file deduped_unique.sam does not match expected output" && exit 1)

############################################################################################

rm -rf "$out_dir"

echo "All tests succeeded!"
exit 0