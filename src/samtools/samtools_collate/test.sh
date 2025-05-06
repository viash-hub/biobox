#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/out"

############################################################################################

echo ">>> Test 1: $meta_name"
"$meta_executable" \
  --input "$test_dir/test.paired_end.sorted.bam" \
  --output "$out_dir/collated.bam"

echo ">>> Checking whether output exists"
[ ! -f "$out_dir/collated.bam" ] && echo "File 'collated.bam' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$out_dir/collated.bam" ] && echo "File 'collated.bam' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff <(samtools view "$out_dir/collated.bam") \
    <(samtools view "$test_dir/collated.bam") || \
    (echo "Output file collated.bam does not match expected output" && exit 1)

############################################################################################

echo ">>> Test 2: $meta_name with --fast option"
"$meta_executable" \
  --fast \
  --input "$test_dir/test.paired_end.sorted.bam" \
  --output "$out_dir/fast_collated.bam"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/fast_collated.bam" ] && echo "File 'fast_collated.bam' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/fast_collated.bam" ] && echo "File 'fast_collated.bam' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff <(samtools view "$test_dir/fast_collated.bam") \
    <(samtools view "$test_dir/fast_collated.bam") || \
    (echo "Output file fast_collated.bam does not match expected output" && exit 1)


############################################################################################

echo ">>> Test 3: $meta_name with compression"
"$meta_executable" \
  --compression 8 \
  --input "$test_dir/test.paired_end.sorted.bam" \
  --output "$out_dir/comp_collated.bam" 

echo ">>> Checking whether output exists"
[ ! -f "$out_dir/comp_collated.bam" ] && echo "File 'comp_collated.bam' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$out_dir/comp_collated.bam" ] && echo "File 'comp_collated.bam' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff <(samtools view "$out_dir/comp_collated.bam") \
    <(samtools view "$test_dir/comp_collated.bam") || \
    (echo "Output file comp_collated.bam does not match expected output" && exit 1)

############################################################################################

echo ">>> All tests passed successfully."

exit 0
