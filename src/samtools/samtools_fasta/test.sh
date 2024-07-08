#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/out_data"

############################################################################################

echo ">>> Test 1: Convert all reads from a bam file to fasta format"
"$meta_executable" \
  --input "$test_dir/a.bam" \
  --output "$out_dir/a.fa"

echo ">>> Check if output file exists"
[ ! -f "$out_dir/a.fa" ] && echo "Output file a.fa does not exist" && exit 1

echo ">>> Check if output is empty"
[ ! -s "$out_dir/a.fa" ] && echo "Output file a.fa is empty" && exit 1

echo ">>> Check if output matches expected output"
diff "$out_dir/a.fa" "$test_dir/a.fa" || 
  (echo "Output file a.fa does not match expected output" && exit 1)

rm "$out_dir/a.fa"

############################################################################################

echo ">>> Test 2: Convert all reads from a sam file to fasta format"
"$meta_executable" \
  --input "$test_dir/a.sam" \
  --output "$out_dir/a.fa"

echo ">>> Check if output file exists"
[ ! -f "$out_dir/a.fa" ] && echo "Output file a.fa does not exist" && exit 1

echo ">>> Check if output is empty"
[ ! -s "$out_dir/a.fa" ] && echo "Output file a.fa is empty" && exit 1

echo ">>> Check if output matches expected output"
diff "$out_dir/a.fa" "$test_dir/a.fa" || 
  (echo "Output file a.fa does not match expected output" && exit 1)

rm "$out_dir/a.fa"

############################################################################################

echo ">>> Test 3: Output reads from bam file to separate files"

"$meta_executable" \
  --input "$test_dir/a.bam" \
  --read1 "$out_dir/a.1.fa" \
  --read2 "$out_dir/a.2.fa" \
  --output "$out_dir/a.fa"

echo ">>> Check if output files exist"
[ ! -f "$out_dir/a.1.fa" ] && echo "Output file a.1.fa does not exist" && exit 1
[ ! -f "$out_dir/a.2.fa" ] && echo "Output file a.2.fa does not exist" && exit 1
[ ! -f "$out_dir/a.fa" ] && echo "Output file a.fa does not exist" && exit 1

echo ">>> Check if output files are empty"
[ ! -s "$out_dir/a.1.fa" ] && echo "Output file a.1.fa is empty" && exit 1
[ ! -s "$out_dir/a.2.fa" ] && echo "Output file a.2.fa is empty" && exit 1
# output should be empty since input has no singleton reads

echo ">>> Check if output files match expected output"
diff "$out_dir/a.1.fa" "$test_dir/a.1.fa" || 
  (echo "Output file a.1.fa does not match expected output" && exit 1)
diff "$out_dir/a.2.fa" "$test_dir/a.2.fa" ||
  (echo "Output file a.2.fa does not match expected output" && exit 1)

rm "$out_dir/a.1.fa" "$out_dir/a.2.fa" "$out_dir/a.fa"

############################################################################################

echo ">>> Test 4: Output only forward reads from bam file to fasta format"

"$meta_executable" \
  --input "$test_dir/a.sam" \
  --excl_flags "0x80" \
  --output "$out_dir/half.fa"

echo ">>> Check if output file exists"
[ ! -f "$out_dir/half.fa" ] && echo "Output file half.fa does not exist" && exit 1

echo ">>> Check if output is empty"
[ ! -s "$out_dir/half.fa" ] && echo "Output file half.fa is empty" && exit 1

echo ">>> Check if output matches expected output"
diff "$out_dir/half.fa" "$test_dir/half.fa" || 
  (echo "Output file half.fa does not match expected output" && exit 1)

rm "$out_dir/half.fa"

############################################################################################

echo "All tests succeeded!"
exit 0