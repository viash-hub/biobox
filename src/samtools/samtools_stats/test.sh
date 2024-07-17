#!/bin/bash

test_dir="${meta_resources_dir}/test_data"

############################################################################################

echo ">>> Test 1: $meta_functionality_name"
"$meta_executable" \
  --input "$test_dir/test.paired_end.sorted.bam" \
  --bai "$test_dir/test.paired_end.sorted.bam.bai" \
  --output "$test_dir/test.paired_end.sorted.txt"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/test.paired_end.sorted.txt" ] && echo "File 'test.paired_end.sorted.txt' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/test.paired_end.sorted.txt" ] && echo "File 'test.paired_end.sorted.txt' is empty!" && exit 1

echo ">>> Checking whether output is correct"
# compare using diff, ignoring the line stating the command that was passed.
diff <(grep -v "^# The command" "$test_dir/test.paired_end.sorted.txt") \
    <(grep -v "^# The command" "$test_dir/ref.paired_end.sorted.txt") || \
    (echo "Output file ref.paired_end.sorted.txt does not match expected output" && exit 1)

rm "$test_dir/test.paired_end.sorted.txt"

############################################################################################

echo ">>> Test 2: $meta_functionality_name with --remove_dups"
"$meta_executable" \
  --remove_dups \
  --input "$test_dir/test.paired_end.sorted.bam" \
  --bai "$test_dir/test.paired_end.sorted.bam.bai" \
  --output "$test_dir/test.d.paired_end.sorted.txt"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/ref.d.paired_end.sorted.txt" ] && echo "File 'ref.d.paired_end.sorted.txt' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/ref.d.paired_end.sorted.txt" ] && echo "File 'ref.d.paired_end.sorted.txt' is empty!" && exit 1

echo ">>> Checking whether output is correct"
# compare using diff, ignoring the line stating the command that was passed.
diff <(grep -v "^# The command" "$test_dir/test.d.paired_end.sorted.txt") \
    <(grep -v "^# The command" "$test_dir/ref.d.paired_end.sorted.txt") || \
    (echo "Output file ref.d.paired_end.sorted.txt does not match expected output" && exit 1)

rm "$test_dir/test.d.paired_end.sorted.txt"

############################################################################################

echo ">>> Test 3: $meta_functionality_name with --remove_overlaps"
"$meta_executable" \
  --remove_overlaps \
  --input "$test_dir/test.paired_end.sorted.bam" \
  --bai "$test_dir/test.paired_end.sorted.bam.bai" \
  --output "$test_dir/test.p.paired_end.sorted.txt" 

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/ref.p.paired_end.sorted.txt" ] && echo "File 'ref.p.paired_end.sorted.txt' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/ref.p.paired_end.sorted.txt" ] && echo "File 'ref.p.paired_end.sorted.txt' is empty!" && exit 1


echo ">>> Checking whether output is correct"
# compare using diff, ignoring the line stating the command that was passed.
diff <(grep -v "^# The command" "$test_dir/test.p.paired_end.sorted.txt") \
    <(grep -v "^# The command" "$test_dir/ref.p.paired_end.sorted.txt") || \
    (echo "Output file ref.p.paired_end.sorted.txt does not match expected output" && exit 1)

rm "$test_dir/test.p.paired_end.sorted.txt"

############################################################################################

echo ">>> All tests passed successfully."

exit 0
