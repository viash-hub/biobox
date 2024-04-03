#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --bam "$test_dir/a.sorted.bam" \
  --bai "$test_dir/a.sorted.bam.bai" \
  --output "$test_dir/a.sorted.idxstats"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/a.sorted.idxstats" ] && echo "File 'a.sorted.idxstats' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/a.sorted.idxstats" ] && echo "File 'a.sorted.idxstats' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$test_dir/a.sorted.idxstats" "$test_dir/a_ref.sorted.idxstats" || \
    (echo "Output file a.sorted.idxstats does not match expected output" && exit 1)

rm "$test_dir/a.sorted.idxstats"

############################################################################################

echo ">>> Testing $meta_functionality_name with singletons in the input"

"$meta_executable" \
  --bam "$test_dir/test.paired_end.sorted.bam" \
  --bai "$test_dir/test.paired_end.sorted.bam.bai" \
  --output "$test_dir/test.paired_end.sorted.idxstats"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/test.paired_end.sorted.idxstats" ] && \
    echo "File 'test.paired_end.sorted.idxstats' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/test.paired_end.sorted.idxstats" ] && \
    echo "File 'test.paired_end.sorted.idxstats' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$test_dir/test.paired_end.sorted.idxstats" "$test_dir/test_ref.paired_end.sorted.idxstats" || \
    (echo "Output file test.paired_end.sorted.idxstats does not match expected output" && exit 1)

rm "$test_dir/test.paired_end.sorted.idxstats"

############################################################################################

echo "All tests succeeded!"
exit 0