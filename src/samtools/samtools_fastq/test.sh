#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/out_data"

############################################################################################

echo ">>> Test 1: Convert all reads from a bam file to fastq format"
"$meta_executable" \
  --input "$test_dir/a.bam" \
  --output "$out_dir/a.fq"

echo ">>> Check if output file exists"
[ ! -f "$out_dir/a.fq" ] && echo "Output file a.fq does not exist" && exit 1

echo ">>> Check if output is empty"
[ ! -s "$out_dir/a.fq" ] && echo "Output file a.fq is empty" && exit 1

echo ">>> Check if output matches expected output"
diff "$out_dir/a.fq" "$test_dir/a.fq" || 
  (echo "Output file a.fq does not match expected output" && exit 1)

rm "$out_dir/a.fq"

############################################################################################

echo ">>> Test 2: Convert all reads from a sam file to fastq format"
"$meta_executable" \
  --input "$test_dir/a.sam" \
  --output "$out_dir/a.fq"

echo ">>> Check if output file exists"
[ ! -f "$out_dir/a.fq" ] && echo "Output file a.fq does not exist" && exit 1

echo ">>> Check if output is empty"
[ ! -s "$out_dir/a.fq" ] && echo "Output file a.fq is empty" && exit 1

echo ">>> Check if output matches expected output"
diff "$out_dir/a.fq" "$test_dir/a.fq" || 
  (echo "Output file a.fq does not match expected output" && exit 1)

rm "$out_dir/a.fq"

############################################################################################

echo ">>> Test 3: Output reads from bam file to separate files"

"$meta_executable" \
  --input "$test_dir/a.bam" \
  --read1 "$out_dir/a.1.fq" \
  --read2 "$out_dir/a.2.fq" \
  --output "$out_dir/a.fq"

echo ">>> Check if output files exist"
[ ! -f "$out_dir/a.1.fq" ] && echo "Output file a.1.fq does not exist" && exit 1
[ ! -f "$out_dir/a.2.fq" ] && echo "Output file a.2.fq does not exist" && exit 1
[ ! -f "$out_dir/a.fq" ] && echo "Output file a.fq does not exist" && exit 1

echo ">>> Check if output files are empty"
[ ! -s "$out_dir/a.1.fq" ] && echo "Output file a.1.fq is empty" && exit 1
[ ! -s "$out_dir/a.2.fq" ] && echo "Output file a.2.fq is empty" && exit 1
# output should be empty since input has no singleton reads

echo ">>> Check if output files match expected output"
diff "$out_dir/a.1.fq" "$test_dir/a.1.fq" || 
  (echo "Output file a.1.fq does not match expected output" && exit 1)
diff "$out_dir/a.2.fq" "$test_dir/a.2.fq" ||
  (echo "Output file a.2.fq does not match expected output" && exit 1)

rm "$out_dir/a.1.fq" "$out_dir/a.2.fq" "$out_dir/a.fq"

############################################################################################

echo ">>> Test 4: Output only forward reads from bam file to fastq format"

"$meta_executable" \
  --input "$test_dir/a.sam" \
  --excl_flags "0x80" \
  --output "$out_dir/half.fq"

echo ">>> Check if output file exists"
[ ! -f "$out_dir/half.fq" ] && echo "Output file half.fq does not exist" && exit 1

echo ">>> Check if output is empty"
[ ! -s "$out_dir/half.fq" ] && echo "Output file half.fq is empty" && exit 1

echo ">>> Check if output matches expected output"
diff "$out_dir/half.fq" "$test_dir/half.fq" || 
  (echo "Output file half.fq does not match expected output" && exit 1)

rm "$out_dir/half.fq"

############################################################################################

echo "All tests succeeded!"
exit 0