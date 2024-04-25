#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/test_data/text"

# Files are compared using the "samtools view" output.
############################################################################################

echo ">>> Test 1: Sorting a BAM file"

"$meta_executable" \
  --input "$test_dir/a.bam" \
  --output "$test_dir/a.sorted.bam"

echo ">>> Check if output file exists"
[ ! -f "$test_dir/a.sorted.bam" ] \
    && echo "Output file a.sorted.bam does not exist" && exit 1

echo ">>> Check if output is empty"
[ ! -s "$test_dir/a.sorted.bam" ] \
    && echo "Output file a.sorted.bam is empty" && exit 1

echo ">>> Check if output matches expected output"
diff -a "$test_dir/a.sorted.bam.txt" "$out_dir/a_ref.sorted.txt" \
    || (echo "Output file a.sorted.bam does not match expected output" && exit 1)

rm "$test_dir/a.sorted.bam" "$test_dir/a.sorted.bam.txt"

############################################################################################

echo ">>> Test 2: Sorting a BAM file according to ascii order"

"$meta_executable" \
  --input "$test_dir/a.bam" \
  --ascii_sort \
  --output "$test_dir/ascii.sorted.bam"

echo ">>> Check if output file exists"
[ ! -f "$test_dir/ascii.sorted.bam" ] \
    && echo "Output file ascii.sorted.bam does not exist" && exit 1

echo ">>> Check if output is empty"
[ ! -s "$test_dir/ascii.sorted.bam" ] \
    && echo "Output file ascii.sorted.bam is empty" && exit 1

echo ">>> Check if output matches expected output"
diff -a "$test_dir/ascii.sorted.bam.txt" "$out_dir/ascii_ref.sorted.txt" \
    || (echo "Output file ascii.sorted.bam does not match expected output" && exit 1)

rm "$test_dir/ascii.sorted.bam" "$test_dir/ascii.sorted.bam.txt"

############################################################################################

echo ">>> Test 3: Sorting a BAM file with compression"

"$meta_executable" \
  --input "$test_dir/a.bam" \
  --compression 5 \
  --output "$test_dir/compressed.sorted.bam"
  
echo ">>> Check if output file exists"
[ ! -f "$test_dir/compressed.sorted.bam" ] \
    && echo "Output file compressed.sorted.bam does not exist" && exit 1

echo ">>> Check if output is empty"
[ ! -s "$test_dir/compressed.sorted.bam" ] \
    && echo "Output file compressed.sorted.bam is empty" && exit 1

echo ">>> Check if output matches expected output" #
diff "$test_dir/compressed.sorted.bam.txt" "$out_dir/compressed_ref.sorted.txt" \
    || (echo "Output file compressed.sorted.bam does not match expected output" && exit 1)

rm "$test_dir/compressed.sorted.bam" "$test_dir/compressed.sorted.bam.txt"

############################################################################################


echo "All tests succeeded!"
exit 0